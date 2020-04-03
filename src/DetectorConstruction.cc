//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes, nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: DetectorConstruction.cc, v 1.18 2010-10-23 19:27:38 gum Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

#include "DetectorConstruction.hh"
#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ClassicalRK4.hh"
#include "G4ExplicitEuler.hh"
#include "G4ChordFinder.hh"
#include "G4EqMagElectricField.hh"
#include "G4PropagatorInField.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SubtractionSolid.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SDManager.hh"

#include "G4UserLimits.hh"

#include "DetectorConstruction.hh"
#include <G4TransportationManager.hh>
#include <G4MagneticField.hh>
#include <G4UniformMagField.hh>
#include <G4FieldManager.hh>
#include "CreateTree.hh"
#include <algorithm>
#include <string>
#include <sstream>
#include "G4MultiUnion.hh"

using namespace CLHEP;

DetectorConstruction::DetectorConstruction(const string &configFileName)
{
  //---------------------------------------
  //------------- Parameters --------------
  //---------------------------------------

  ConfigFile config(configFileName);

  config.readInto(checkOverlaps, "checkOverlaps");

  config.readInto(world_material, "world_material");
  config.readInto(bar_length, "bar_length");

  config.readInto(core_radius_x, "core_radius_x");
  config.readInto(core_radius_y, "core_radius_y");
  config.readInto(core_material, "core_material");
  config.readInto(core_rIndex, "core_rIndex");
  config.readInto(core_absLength, "core_absLength");

  config.readInto(gap_l, "gap_l");
  config.readInto(gap_size_x, "gap_size_x");
  config.readInto(gap_size_y, "gap_size_y");
  config.readInto(gap_material, "gap_material");

  config.readInto(det_l, "det_l");
  config.readInto(det_size_x, "det_size_x");
  config.readInto(det_size_y, "det_size_y");
  config.readInto(det_material, "det_material");

  config.readInto(depth, "depth");
  config.readInto(cryst_dist, "cryst_dist");
  config.readInto(trackerX0, "trackerX0");
  config.readInto(services_thick, "services_thick");

  config.readInto(ecal_incline, "ecal_incline");
  config.readInto(ecal_material, "ecal_material");
  config.readInto(scinti_material, "scinti_material");
  config.readInto(Cherenc_material, "Cherenc_material");
  config.readInto(Cherenp_material, "Cherenp_material");

  config.readInto(ecal_front_length, "ecal_front_length");
  config.readInto(ecal_rear_length, "ecal_rear_length");
  config.readInto(ecal_front_face, "ecal_front_face");
  config.readInto(ecal_rear_face, "ecal_rear_face");
  config.readInto(hole_diameter, "hole_diameter");
  config.readInto(fiber_diameter, "fiber_diameter");
  config.readInto(ecal_timing_distance, "ecal_timing_distance");
  config.readInto(ecal_det_size, "ecal_det_size");

  config.readInto(hcal_width, "hcal_width");
  config.readInto(hcalTile_width, "hcalTile_width");
  config.readInto(hcalAbs_1_thick, "hcalAbs_1_thick");
  config.readInto(hcalAbs_2_thick, "hcalAbs_2_thick");
  config.readInto(solenoid_thick, "solenoid_thick");
  config.readInto(hcalTile_thick, "hcalTile_thick");

  B_field_intensity = config.read<double>("B_field_intensity") * tesla;

  expHall_x = 300. * cm;
  expHall_y = 300. * cm;
  expHall_z = 1500. * cm;

  B_field_IsInitialized = false;

  initializeMaterials();

  CreateTree::Instance()->inputTrackerX0 = trackerX0;
  CreateTree::Instance()->inputServiceAlmm = services_thick;
  CreateTree::Instance()->inputTimingThick = core_radius_x * 2;
  CreateTree::Instance()->inputE1Thick = ecal_front_length;
  CreateTree::Instance()->inputE2Thick = ecal_rear_length;
  CreateTree::Instance()->inputE1Width = ecal_front_face;
  CreateTree::Instance()->inputTimingECAL_dist = ecal_timing_distance;
}

//---- ---- ---- ---- ---- ---- ---- ---- ----  ---- ---- ---- ---- ---- ----

DetectorConstruction::~DetectorConstruction()
{
  delete stepLimit;
}

//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

G4VPhysicalVolume *DetectorConstruction::Construct()
{
  G4cout << ">>>>>> DetectorConstruction::Construct ()::begin <<<<<<" << G4endl;

  //------------------------------------
  //------------- Geometry -------------
  //------------------------------------

  // The experimental Hall
  // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  G4VSolid *worldS = new G4Box("worldS", 5 * expHall_x, 5 * expHall_y, 5 * expHall_z);
  G4LogicalVolume *worldLV = new G4LogicalVolume(worldS, WoMaterial, "worldLV", 0, 0, 0);
  G4VPhysicalVolume *worldPV = new G4PVPlacement(0, G4ThreeVector(), worldLV, "worldPV", 0, false, 0, checkOverlaps);


  // the cooling-services Dead material layer: 5 cm
  // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


  const int NECAL_CRYST = 100; 




  //********************************************
  //  ELECTROMAGNETIC CALORIMETER
  //********************************************
  std::cout<<"EM calorimeter construction "<<std::endl;



  std::cout<<"making passive with width "<<hole_diameter<<std::endl;
  G4Box *ecalCrystalS_f_absorb = new G4Box("ecalCrystalS_f_absorb", 0.5*ecal_front_face, 0.5 * ecal_rear_face,  0.5*(hole_diameter-5.));

  G4LogicalVolume *ecalCrystalL_f_absorb = new G4LogicalVolume(ecalCrystalS_f_absorb, EcalMaterial, "ecalCrystalL_f_absorb", 0, 0, 0);



  std::cout<<"making active with width "<<fiber_diameter<<std::endl;
  G4Box *ecalCrystalS_f_fiber = new G4Box("ecalCrystalS_f_fiber", 0.5*ecal_front_face, 0.5 * ecal_rear_face, 0.5 * (3.*fiber_diameter));

  G4LogicalVolume *ecalCrystalL_f_fiber_scinti = new G4LogicalVolume(ecalCrystalS_f_fiber, ScintiMaterial, "ecalCrystalL_f_fiber_scinti", 0, 0, 0);
  //  G4LogicalVolume *ecalCrystalL_f_fiber_cherenc = new G4LogicalVolume(ecalCrystalS_f_fiber, CherencMaterial, "ecalCrystalL_f_fiber_cherenc", 0, 0, 0);
  //  G4LogicalVolume *ecalCrystalL_f_fiber_cherenp = new G4LogicalVolume(ecalCrystalS_f_fiber, CherenpMaterial, "ecalCrystalL_f_fiber_cherenp", 0, 0, 0);


  //


  // ECAL physical placement
  G4VPhysicalVolume *ecalCrystalP_f[NECAL_CRYST];




  G4VPhysicalVolume *ecalCrystalP_f_fiber_scinti[NECAL_CRYST];


  //G4VPhysicalVolume *ecalCrystalP_f_fiber_cherenc[NECAL_CRYST];


  //G4VPhysicalVolume *ecalCrystalP_f_fiber_cherenp[NECAL_CRYST];




  char name[60];

  G4double z_pos_f[NECAL_CRYST];

  int nArrayECAL = (int)sqrt(NECAL_CRYST);

  int iCrystal;
  for (int iX = 0; iX < nArrayECAL; iX++)
  {
    for (int iY = 0; iY < nArrayECAL; iY++)
    {

      G4RotationMatrix *piRotEcal = new G4RotationMatrix;
      //piRotEcal->rotateX(pointingAngle * deg);

      iCrystal = nArrayECAL * iX + iY;

      z_pos_f[iCrystal] = iCrystal*( hole_diameter+ fiber_diameter*3);

      std::cout<<"iCrystal "<<iCrystal<<" iX iY "<<iX<<" "<<iY<<" z "<<z_pos_f[iCrystal]<<std::endl;



      //add dream detector instead of the original one

      sprintf(name, "ecalCrystalP_f_absorb_%d", iCrystal);
      ecalCrystalP_f[iCrystal] = new G4PVPlacement(0, G4ThreeVector(0, 0, ecal_timing_distance + z_pos_f[iCrystal]+0.5*hole_diameter), ecalCrystalL_f_absorb, name, worldLV, false, 0);

      sprintf(name, "ecalCrystalP_f_fiber_scinti_%d", iCrystal);
      ecalCrystalP_f_fiber_scinti[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(0, 0, ecal_timing_distance + z_pos_f[iCrystal]+hole_diameter+0.5*3.*fiber_diameter), ecalCrystalL_f_fiber_scinti, name, worldLV, false, 0);

      /*
      sprintf(name, "ecalCrystalP_f_fiber_cherenc_%d", iCrystal);
      ecalCrystalP_f_fiber_cherenc[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(0, 0, ecal_timing_distance + z_pos_f[iCrystal]+hole_diameter+1.5*fiber_diameter), ecalCrystalL_f_fiber_cherenc, name, worldLV, false, 0);


      sprintf(name, "ecalCrystalP_f_fiber_cherenp_%d", iCrystal);
      ecalCrystalP_f_fiber_cherenp[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(0, 0, ecal_timing_distance + z_pos_f[iCrystal]+hole_diameter+2.5*fiber_diameter), ecalCrystalL_f_fiber_cherenp, name, worldLV, false, 0);
      */

      //

      sprintf(name, "ecalGapP_f_%d", iCrystal);




      sprintf(name, "ecalDetP_f_%d", iCrystal);




    }
  }


  //-----------------------------------------------------
  //------------- Visualization attributes --------------
  //-----------------------------------------------------

  G4Colour white(1.00, 1.00, 1.00);   // white
  G4Colour gray(0.50, 0.50, 0.50);    // gray
  G4Colour black(0.00, 0.00, 0.00);   // black
  G4Colour red(1.00, 0.00, 0.00);     // red
  G4Colour green(0.00, 1.00, 0.00);   // green
  G4Colour blue(0.00, 0.00, 1.00);    // blue
  G4Colour cyan(0.00, 1.00, 1.00);    // cyan
  G4Colour air(0.90, 0.94, 1.00);     // cyan
  G4Colour magenta(1.00, 0.00, 1.00); // magenta
  G4Colour yellow(1.00, 1.00, 0.00);  // yellow
  G4Colour brass(0.80, 0.60, 0.40);   // brass
  G4Colour brown(0.70, 0.40, 0.10);   // brown

  G4VisAttributes *VisAttWorld = new G4VisAttributes(black);
  VisAttWorld->SetVisibility(true);
  VisAttWorld->SetForceWireframe(true);
  worldLV->SetVisAttributes(VisAttWorld);






  G4VisAttributes *VisCrystalCore = new G4VisAttributes(blue);
  VisCrystalCore->SetVisibility(true);
  VisCrystalCore->SetForceWireframe(true);

  ecalCrystalL_f_absorb->SetVisAttributes(VisCrystalCore);

  G4VisAttributes *VisFiber = new G4VisAttributes(red);
  VisFiber->SetVisibility(true);
  VisFiber->SetForceWireframe(true);
  ecalCrystalL_f_fiber_scinti->SetVisAttributes(VisFiber);
  //  ecalCrystalL_f_fiber_cherenc->SetVisAttributes(VisFiber);
  //ecalCrystalL_f_fiber_cherenp->SetVisAttributes(VisFiber);


  if (B_field_intensity > 0.1 * tesla)
    ConstructField();


  G4cout << ">>>>>> DetectorConstruction::Construct ()::end <<< " << G4endl;
  return worldPV;
}

//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

void DetectorConstruction::initializeMaterials()
{
  //-----------------
  // define materials

  std::cout<<"entering initializeMaterials"<<std::endl;

  WoMaterial = NULL;
  if (world_material == 0)
    WoMaterial = MyMaterials::Vacuum();
  else if (world_material == 1)
    WoMaterial = MyMaterials::Air();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre world material specifier " << world_material << G4endl;
    exit(-1);
  }
  G4cout << "Wo. material: " << WoMaterial << G4endl;

  CoMaterial = NULL;
  if (core_material == 1)
    CoMaterial = MyMaterials::Quartz();
  else if (core_material == 2)
    CoMaterial = MyMaterials::SiO2();
  else if (core_material == 3)
    CoMaterial = MyMaterials::SiO2_Ce();
  else if (core_material == 4)
    CoMaterial = MyMaterials::LuAG_Ce();
  else if (core_material == 5)
    CoMaterial = MyMaterials::YAG_Ce();
  else if (core_material == 6)
    CoMaterial = MyMaterials::LSO();
  else if (core_material == 7)
    CoMaterial = MyMaterials::LYSO();
  else if (core_material == 8)
    CoMaterial = MyMaterials::LuAG_undoped();
  else if (core_material == 9)
    CoMaterial = MyMaterials::GAGG_Ce();
  else if (core_material == 11)
    CoMaterial = MyMaterials::LuAG_Pr();
  else if (core_material == 12)
    CoMaterial = MyMaterials::PbF2();
  else if (core_material == 13)
    CoMaterial = MyMaterials::PlasticBC408();
  else if (core_material == 14)
    CoMaterial = MyMaterials::PlasticBC418();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << core_material << G4endl;
    exit(-1);
  }
  G4cout << "Co. material: " << CoMaterial << G4endl;

  EcalMaterial = NULL;
  if (ecal_material == 1)
    EcalMaterial = MyMaterials::Quartz();
  else if (ecal_material == 2)
    EcalMaterial = MyMaterials::SiO2();
  else if (ecal_material == 3)
    EcalMaterial = MyMaterials::SiO2_Ce();
  else if (ecal_material == 4)
    EcalMaterial = MyMaterials::LuAG_Ce();
  else if (ecal_material == 5)
    EcalMaterial = MyMaterials::YAG_Ce();
  else if (ecal_material == 6)
    EcalMaterial = MyMaterials::LSO();
  else if (ecal_material == 7)
    EcalMaterial = MyMaterials::LYSO();
  else if (ecal_material == 8)
    EcalMaterial = MyMaterials::LuAG_undoped();
  else if (ecal_material == 9)
    EcalMaterial = MyMaterials::GAGG_Ce();
  else if (ecal_material == 10)
    EcalMaterial = MyMaterials::LuAG_Pr();
  else if (ecal_material == 11)
    EcalMaterial = MyMaterials::PbF2();
  else if (ecal_material == 12)
    EcalMaterial = MyMaterials::PlasticBC408();
  else if (ecal_material == 13)
    EcalMaterial = MyMaterials::PlasticBC418();
  else if (ecal_material == 14)
    EcalMaterial = MyMaterials::PWO();
  else if (ecal_material == 15)
    EcalMaterial = MyMaterials::Acrylic();
  else if (ecal_material == 16)
    EcalMaterial = MyMaterials::copper();
  else if (ecal_material == 17)
    EcalMaterial = MyMaterials::Brass();
  else if (ecal_material == 18)
    EcalMaterial = MyMaterials::EJ200();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << ecal_material << G4endl;
    exit(-1);
  }
  G4cout << "ECAL material: " << EcalMaterial << G4endl;

  /************************************************************************************/
  ScintiMaterial = NULL;
  if (scinti_material == 1)
    ScintiMaterial = MyMaterials::Quartz();
  else if (scinti_material == 2)
    ScintiMaterial = MyMaterials::SiO2();
  else if (scinti_material == 3)
    ScintiMaterial = MyMaterials::SiO2_Ce();
  else if (scinti_material == 4)
    ScintiMaterial = MyMaterials::LuAG_Ce();
  else if (scinti_material == 5)
    ScintiMaterial = MyMaterials::YAG_Ce();
  else if (scinti_material == 6)
    ScintiMaterial = MyMaterials::LSO();
  else if (scinti_material == 7)
    ScintiMaterial = MyMaterials::LYSO();
  else if (scinti_material == 8)
    ScintiMaterial = MyMaterials::LuAG_undoped();
  else if (scinti_material == 9)
    ScintiMaterial = MyMaterials::GAGG_Ce();
  else if (scinti_material == 10)
    ScintiMaterial = MyMaterials::LuAG_Pr();
  else if (scinti_material == 11)
    ScintiMaterial = MyMaterials::PbF2();
  else if (scinti_material == 12)
    ScintiMaterial = MyMaterials::PlasticBC408();
  else if (scinti_material == 13)
    ScintiMaterial = MyMaterials::PlasticBC418();
  else if (scinti_material == 14)
    ScintiMaterial = MyMaterials::PWO();
  else if (scinti_material == 15)
    ScintiMaterial = MyMaterials::Acrylic();
  else if (scinti_material == 16)
    ScintiMaterial = MyMaterials::copper();
  else if (scinti_material == 17)
    ScintiMaterial = MyMaterials::Brass();
  else if (scinti_material == 18)
    ScintiMaterial = MyMaterials::EJ200();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << scinti_material << G4endl;
    exit(-1);
  }
  G4cout << "SCINT material: " << ScintiMaterial << G4endl;


  CherencMaterial = NULL;
  if (Cherenc_material == 1)
    CherencMaterial = MyMaterials::Quartz();
  else if (Cherenc_material == 2)
    CherencMaterial = MyMaterials::SiO2();
  else if (Cherenc_material == 3)
    CherencMaterial = MyMaterials::SiO2_Ce();
  else if (Cherenc_material == 4)
    CherencMaterial = MyMaterials::LuAG_Ce();
  else if (Cherenc_material == 5)
    CherencMaterial = MyMaterials::YAG_Ce();
  else if (Cherenc_material == 6)
    CherencMaterial = MyMaterials::LSO();
  else if (Cherenc_material == 7)
    CherencMaterial = MyMaterials::LYSO();
  else if (Cherenc_material == 8)
    CherencMaterial = MyMaterials::LuAG_undoped();
  else if (Cherenc_material == 9)
    CherencMaterial = MyMaterials::GAGG_Ce();
  else if (Cherenc_material == 10)
    CherencMaterial = MyMaterials::LuAG_Pr();
  else if (Cherenc_material == 11)
    CherencMaterial = MyMaterials::PbF2();
  else if (Cherenc_material == 12)
    CherencMaterial = MyMaterials::PlasticBC408();
  else if (Cherenc_material == 13)
    CherencMaterial = MyMaterials::PlasticBC418();
  else if (Cherenc_material == 14)
    CherencMaterial = MyMaterials::PWO();
  else if (Cherenc_material == 15)
    CherencMaterial = MyMaterials::Acrylic();
  else if (Cherenc_material == 16)
    CherencMaterial = MyMaterials::copper();
  else if (Cherenc_material == 17)
    CherencMaterial = MyMaterials::Brass();
  else if (Cherenc_material == 18)
    CherencMaterial = MyMaterials::EJ200();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << Cherenc_material << G4endl;
    exit(-1);
  }
  G4cout << "Cherenk material: " << CherencMaterial << G4endl;


  CherenpMaterial = NULL;
  if (Cherenp_material == 1)
    CherenpMaterial = MyMaterials::Quartz();
  else if (Cherenp_material == 2)
    CherenpMaterial = MyMaterials::SiO2();
  else if (Cherenp_material == 3)
    CherenpMaterial = MyMaterials::SiO2_Ce();
  else if (Cherenp_material == 4)
    CherenpMaterial = MyMaterials::LuAG_Ce();
  else if (Cherenp_material == 5)
    CherenpMaterial = MyMaterials::YAG_Ce();
  else if (Cherenp_material == 6)
    CherenpMaterial = MyMaterials::LSO();
  else if (Cherenp_material == 7)
    CherenpMaterial = MyMaterials::LYSO();
  else if (Cherenp_material == 8)
    CherenpMaterial = MyMaterials::LuAG_undoped();
  else if (Cherenp_material == 9)
    CherenpMaterial = MyMaterials::GAGG_Ce();
  else if (Cherenp_material == 10)
    CherenpMaterial = MyMaterials::LuAG_Pr();
  else if (Cherenp_material == 11)
    CherenpMaterial = MyMaterials::PbF2();
  else if (Cherenp_material == 12)
    CherenpMaterial = MyMaterials::PlasticBC408();
  else if (Cherenp_material == 13)
    CherenpMaterial = MyMaterials::PlasticBC418();
  else if (Cherenp_material == 14)
    CherenpMaterial = MyMaterials::PWO();
  else if (Cherenp_material == 15)
    CherenpMaterial = MyMaterials::Acrylic();
  else if (Cherenp_material == 16)
    CherenpMaterial = MyMaterials::copper();
  else if (Cherenp_material == 17)
    CherenpMaterial = MyMaterials::Brass();
  else if (Cherenp_material == 18)
    CherenpMaterial = MyMaterials::EJ200();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << Cherenp_material << G4endl;
    exit(-1);
  }
  G4cout << "Cherenk2 material: " << CherenpMaterial << G4endl;
  /************************************************************************************/

  GaMaterial = NULL;
  if (gap_material == 1)
    GaMaterial = MyMaterials::Air();
  else if (gap_material == 2)
    GaMaterial = MyMaterials::OpticalGrease();
  else if (gap_material == 3)
    GaMaterial = MyMaterials::MeltMount168();
  else if (gap_material == 4)
    GaMaterial = MyMaterials::OpticalGrease155();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid gap material specifier " << gap_material << G4endl;
    exit(-1);
  }
  G4cout << "Gap material: " << GaMaterial << G4endl;

  DeMaterial = NULL;
  if (det_material == 1)
    DeMaterial = MyMaterials::Silicon();
  else if (det_material == 2)
    DeMaterial = MyMaterials::Quartz();
  else if (det_material == 3)
    DeMaterial = MyMaterials::Air();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid detector material specifier " << det_material << G4endl;
    exit(-1);
  }
  G4cout << "Detector material: " << DeMaterial << G4endl;

  //------------------
  // change properties

  if (core_absLength > 0)
  {
    const G4int nEntries_ABS = 2;
    G4double PhotonEnergy_ABS[nEntries_ABS] = {1. * eV, 10. * eV};
    G4double Absorption[nEntries_ABS] = {core_absLength * mm, core_absLength * mm};

    CoMaterial->GetMaterialPropertiesTable()->RemoveProperty("ABSLENGTH");
    CoMaterial->GetMaterialPropertiesTable()->AddProperty("ABSLENGTH", PhotonEnergy_ABS, Absorption, nEntries_ABS);
  }
  if (core_rIndex > 0)
  {
    const G4int nEntries_RI = 2;
    G4double PhotonEnergy_RI[nEntries_RI] = {1. * eV, 10. * eV};
    G4double RefractiveIndex[nEntries_RI] = {core_rIndex, core_rIndex};

    CoMaterial->GetMaterialPropertiesTable()->RemoveProperty("RINDEX");
    CoMaterial->GetMaterialPropertiesTable()->AddProperty("RINDEX", PhotonEnergy_RI, RefractiveIndex, nEntries_RI);
  }
}

//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

void DetectorConstruction::ConstructField()
{
  G4cout << ">>>>>> DetectorConstruction::ConstructField ()::begin <<<<<<" << G4endl;

  static G4TransportationManager *trMgr = G4TransportationManager::GetTransportationManager();

  // A field object is held by a field manager
  // Find the global Field Manager
  G4FieldManager *globalFieldMgr = trMgr->GetFieldManager();

  if (!B_field_IsInitialized)
  {
    // magnetic field parallel to the beam direction (w/ tilt)
    G4ThreeVector fieldVector(0.0522 * B_field_intensity, 0.0522 * B_field_intensity, 0.9973 * B_field_intensity);

    B_field = new G4UniformMagField(fieldVector);
    globalFieldMgr->SetDetectorField(B_field);
    globalFieldMgr->CreateChordFinder(B_field);
    globalFieldMgr->GetChordFinder()->SetDeltaChord(0.005 * mm);
    B_field_IsInitialized = true;
  }

  G4cout << ">>>>>> DetectorConstruction::ConstructField ()::end <<< " << G4endl;
  return;
}

void DetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((stepLimit) && (maxStep > 0.))
    stepLimit->SetMaxAllowedStep(maxStep);
}