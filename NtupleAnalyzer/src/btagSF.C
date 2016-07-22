#include "BTagCalibrationXStandalone.cpp"
#include "../include/BTagCalibrationXStandalone.h"

void btagSF2(){
    BTagCalibrationX calib = BTagCalibrationX("csvv2",
            "/opt/sbg/scratch1/cms/TTH/weight/CSVv2_ichep.csv");

    BTagCalibrationXReader reader = BTagCalibrationXReader(  BTagEntryX::OP_LOOSE,  // operating point
            "central",            // central sys type
            {"up", "down"});      // other sys types

    reader.load(    calib,                // calibration instance
            BTagEntryX::FLAV_B,    // btag flavour
            "comb");               // measurement type

    std::cout << "So far, so good" << std::endl;

    double jet_scalefactor = 1.;
    jet_scalefactor = reader.eval_auto_bounds( "central", BTagEntryX::FLAV_B, 0.25, 200.2);
    std::cout << jet_scalefactor << std::endl;
}
