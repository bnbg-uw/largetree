// largetree.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "ReadProcessedFolder.hpp"
#include "taolist.hpp"
#include <boost/program_options.hpp>
#include <chrono>

template<class E>
void print(const E& e) {
    std::cout << e << '\n';
}

void silenceGDALErrors(CPLErr eErrClass, CPLErrorNum nError, const char* pszErrorMsg) {
}

int main(int argc, char* argv[])
{
    namespace po = boost::program_options;

    po::options_description desc("LargeTrees: Given lidar data, write out LMU rasters for the lidar data.");
    desc.add_options()
        ("lidar,l", po::value<std::string>(), "REQUIRED: Path to lidar dataset. Expects the FINAL_ folder from fusion. Needs layout files and topometrics.")
        ("fia,f", po::value<std::string>(), "REQUIRED: Path to fia data.")
        ("clim,c", po::value<std::string>(), "REQUIRED: Clim class layers.")
        ("thread,t", po::value<int>(), "how many theads?")
        ("mask,m", po::value<std::string>(), "Mask raster.")
        ("dia,f", po::value<double>(), "diameter cutoff of interest in centimeters. default 76.2 (30inches)")
        ("help,h", "Display this help message and exit.");

    if (argc == 1) {
        std::cout << desc << "\n";
        return 0;
    }

    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
    }
    catch (const po::error& e) {
        std::cerr << "Error: " << e.what() << "\n";
        std::cerr << desc << "\n";
        return 1;
    }

    if (vm.count("help")) {
        std::cout << desc << '\n';
        return 0;
    }

    if (!vm.count("lidar")) {
        std::cerr << "Some required arguments are missing.\n";
        return 1;
    }

    int nthread = 1;
    if (vm.count("thread")) {
        nthread = vm["thread"].as<int>();
    }

    std::unique_ptr<processedfolder::ProcessedFolder> lds;
    lds = processedfolder::readProcessedFolder(vm["lidar"].as<std::string>());
    std::cout << "LDS made...";
    
    auto maskr = lapis::Raster<int>(processedfolder::stringOrThrow(lds->maskRaster()));
    if (vm.count("mask")) {
        auto m = lapis::Raster<int>(vm["mask"].as<std::string>());
        m = lapis::resampleRaster(m, maskr, lapis::ExtractMethod::near);
        m.mask(maskr);
        maskr = m;
    }

    lapis::Raster<int> climate;
    try {
        climate = lapis::Raster<int>(vm["clim"].as<std::string>());
        std::cout << "climate loaded...";
    }
    catch (lapis::InvalidRasterFileException e) {
        std::cerr << "Unable to load climate layer from specified file";
        return 1; 
    }

    auto output = lds->dir();
    if (std::filesystem::exists(output)) {
        output = output / "largetree";
        std::filesystem::create_directory(output);
    } 
    else {
        throw std::runtime_error("Output path does not exist.");
    }

    double dia = 76.2;
    if (vm.count("dia")) {
        dia = vm["dia"].as<double>();
    }
    double fixed = 3;

    std::cout << "Reading TAOs\n";
    auto before = std::chrono::high_resolution_clock::now();
    auto tl = rxtools::TaoListMP(lds->allHighPoints(nthread, fixed),);

    auto after = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(after - before);
    print(duration.count());

    // Done loading input stuff
    //starting fia models

    std::string fiafolder = vm["fia"].as<std::string>();
    auto allPlots = stats::FIALeastSquares::getPlotList(fiafolder);

    lapis::coord_t buffer = 2000 / lds->getConvFactor();
    lapis::Extent bufferext = lapis::Extent(maskr.xmin() - buffer, maskr.xmax() + buffer, maskr.ymin() - buffer, maskr.ymax() + buffer);
    bufferext.crs(maskr.crs());

    stats::FIALeastSquares basedbh = stats::FIALeastSquares(allPlots, fiafolder, "DIA");
    std::size_t plotcount = basedbh.limitByExtent(bufferext);
    if (plotcount == 0) {
        std::cerr << "No FIA plots within vicinity of acquisition. Cannot calculate basal area or biomass.\n";
        exit(EXIT_FAILURE);
    }
    std::cout << "Class-Agnostic:\n";
    std::cout << "DBH:\n";
    basedbh.initializeModel();
    std::cout << basedbh << "\n";
    //climate individual models.
    std::unordered_map<int, stats::FIALeastSquares> dbhByClim;
    std::unordered_set<int> checkedValues;

    for (lapis::cell_t cell = 0; cell < climate.ncell(); ++cell) {
        auto v = climate[cell];
        if (!v.has_value()) {
            continue;
        }
        if (checkedValues.count(v.value())) {
            continue;
        }
        checkedValues.insert(v.value());

        stats::FIALeastSquares thisDbh = basedbh;
        plotcount = thisDbh.limitByRasterValue(climate, v.value());
        if (plotcount > 0) {
            dbhByClim.emplace(v.value(), thisDbh);
            std::cout << "Class " << v.value() << ":\n";
            std::cout << "DBH:\n";
            dbhByClim[v.value()].initializeModel();
            std::cout << dbhByClim[v.value()] << "\n";
        }
    }
    std::cout << "making layers\n";

    lapis::Raster<int> ntree{ lapis::Alignment(maskr) };
    lapis::Raster<int> ntree90{ lapis::Alignment(maskr) };
    lapis::Raster<double> area{ lapis::Alignment(maskr) };
    lapis::Raster<double> area90{ lapis::Alignment(maskr) };
    lapis::Raster<double> totarea{ lapis::Alignment(maskr) };
    lapis::Raster<double> totarea90{ lapis::Alignment(maskr) };
    lapis::Raster<int> ltree{ lapis::Alignment(maskr) };
    lapis::Raster<int> ltree90{ lapis::Alignment(maskr) };
    
    auto fveg = lapis::Raster<int>("F:/fvegyose.img");
    auto cpad = lapis::Raster<int>("F:/cpadyose.img");
    std::ofstream f;
    f.open("F:/trees.csv");
    auto projclim = lapis::resampleRaster(climate, maskr, lapis::ExtractMethod::near);
    for (size_t t = 0; t < tl.size(); ++t) {
        auto x = tl.x()[t];
        auto y = tl.y()[t];
        auto h = tl.height()[t];
        auto thisClim = projclim.extract(x, y, lapis::ExtractMethod::near);
        double dbh = 0;
        if (thisClim.has_value()) {
            if (dbhByClim.count(thisClim.value())) {
                dbh = dbhByClim[thisClim.value()].predictAsMetric(h);
            }
            else {
                dbh = basedbh.predictAsMetric(h);
            }
        }
        else {
            dbh = basedbh.predictAsMetric(h);
        }

        if (!maskr.contains(x, y)) {
            continue;
        }
        if (fveg.extract(x, y, lapis::ExtractMethod::near).has_value() && cpad.extract(x, y, lapis::ExtractMethod::near).has_value()) {
            f << dbh << "," << h << "," << x << "," << y << "\n";
        }
        lapis::rowcol_t row = maskr.rowFromY(y);
        lapis::rowcol_t col = maskr.colFromX(x);
        ntree.atRC(row, col).has_value() = true;
        ntree.atRC(row, col).value()++;
        totarea.atRC(row, col).has_value() = true;
        totarea.atRC(row, col).value() += tl.area()[t];
        ltree.atRC(row, col).has_value() = true;
        if (dbh > dia) {
            ltree.atRC(row, col).value()++;
            area.atRC(row, col).has_value() = true;
            area.atRC(row, col).value() += tl.area()[t];
        }

        for (int rowbudge : {-1, 0, 1}) {
            lapis::rowcol_t thisRow = row + rowbudge;
            if (thisRow < 0 || thisRow >= maskr.nrow()) {
                continue;
            }
            for (int colbudge : {-1, 0, 1}) {
                lapis::rowcol_t thisCol = col + colbudge;
                if (thisCol < 0 || thisCol >= maskr.ncol()) {
                    continue;
                }
                ntree90.atRC(thisRow, thisCol).has_value() = true;
                ntree90.atRC(thisRow, thisCol).value()++;
                totarea90.atRC(row, col).has_value() = true;
                totarea90.atRC(row, col).value() += tl.area()[t];
                ltree90.atRC(thisRow, thisCol).has_value() = true;
                if (dbh > dia) {
                    ltree90.atRC(thisRow, thisCol).value()++;
                    area90.atRC(row, col).has_value() = true;
                    area90.atRC(row, col).value() += tl.area()[t];
                }
            }
        }

    }
    f.close();

    auto tph = ntree / (ntree.xres() * ntree.yres()) * 10000.;
    auto tph90 = ntree90 / (3*ntree.xres() * 3*ntree.yres()) * 10000.;

    std::filesystem::path outfile = output / ("ntree.img");
    ntree.writeRaster(outfile.string());
    outfile = output / ("ntree90.img");
    ntree90.writeRaster(outfile.string());
    outfile = output / ("ltree.img");
    ltree.writeRaster(outfile.string());
    outfile = output / ("ltree90.img");
    ltree90.writeRaster(outfile.string());
    outfile = output / ("tph.img");
    tph.writeRaster(outfile.string());
    outfile = output / ("tph90.img");
    tph90.writeRaster(outfile.string());

    area /= totarea;
    area90 /= totarea90;
    auto all = lapis::Raster<int>(processedfolder::stringOrThrow(lds->maskRaster()));
    for (lapis::cell_t c = 0; c < area.ncell(); ++c) {
        if (!area[c].has_value() && all[c].has_value()) {
            area[c].has_value() = true;
        }
        if (!area90[c].has_value() && all[c].has_value()) {
            area90[c].has_value() = true;
        }
    }

    outfile = output / ("ltreearea.img");
    area.writeRaster(outfile.string());
    outfile = output / ("ltreearea90.img");
    area90.writeRaster(outfile.string());

    auto ltreeph = ltree / (ltree.xres() * ltree.yres()) * 10000.;
    outfile = output / "ltreeph.img";
    ltreeph.writeRaster(outfile.string());

    auto ltreeph90 = ltree90 / (3 * ltree.xres() * 3 * ltree.yres()) * 10000.;
    outfile = output / "ltreeph90.img";
    ltreeph90.writeRaster(outfile.string());

    auto percentltreeph = ltreeph / tph;
    outfile = output / "perccentltreeph.img";
    percentltreeph.writeRaster(outfile.string());

    auto percentltreeph90 = ltreeph90 / tph90;
    outfile = output / "percentltreeph90.img";
    percentltreeph90.writeRaster(outfile.string());

    for (int i = 10; i <= 50; i = i + 10) {
        auto patch = ltreeph90 > i;
        for (lapis::cell_t c = 0; c < patch.ncell(); ++c) {
            patch[c].has_value() = patch[c].value();
        }

        auto cc = lapis::connectedComponents(patch, false);
        std::unordered_map<int, int> regionType;
        for (lapis::cell_t c = 0; c < cc.ncell(); ++c) {
            if (cc[c].has_value()) {
                regionType[cc[c].value()]++;
            }
        }
        
        outfile = output / ("patchid90_" + std::to_string(i) + ".img");
        cc.writeRaster(outfile.string());
    }
    return EXIT_SUCCESS;
}