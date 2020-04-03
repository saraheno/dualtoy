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


  config.readInto(ecal_material, "ecal_material");
  config.readInto(scinti_material, "scinti_material");
  config.readInto(Cherenc_material, "Cherenc_material");


  config.readInto(ecal_front_length, "ecal_front_length");
  config.readInto(ecal_front_face, "ecal_front_face");
  config.readInto(ecal_rear_face, "ecal_rear_face");
  config.readInto(hole_diameter, "hole_diameter");
  config.readInto(fiber_diameter, "fiber_diameter");
  config.readInto(ecal_det_size, "ecal_det_size");



  B_field_intensity = config.read<double>("B_field_intensity") * tesla;

  expHall_x = 300. * cm;
  expHall_y = 300. * cm;
  expHall_z = 1500. * cm;

  B_field_IsInitialized = false;

  initializeMaterials();

  //CreateTree::Instance()->inputTrackerX0 = trackerX0;
  //CreateTree::Instance()->inputServiceAlmm = services_thick;
  //  CreateTree::Instance()->inputTimingThick = core_radius_x * 2;
  CreateTree::Instance()->inputE1Thick = ecal_front_length;
  //CreateTree::Instance()->inputE2Thick = ecal_rear_length;
  CreateTree::Instance()->inputE1Width = ecal_front_face;
  //CreateTree::Instance()->inputTimingECAL_dist = ecal_timing_distance;
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




  const int NECAL_CRYST = 100; 




  //********************************************
  // CALORIMETER
  //********************************************
  std::cout<<"calorimeter construction "<<std::endl;



  std::cout<<"making passive with width "<<hole_diameter<<std::endl;
  G4Box *ecalCrystalS_f_absorb = new G4Box("ecalCrystalS_f_absorb", 0.5*ecal_front_face, 0.5 * ecal_rear_face,  0.5*(hole_diameter-5.));
  G4LogicalVolume *ecalCrystalL_f_absorb = new G4LogicalVolume(ecalCrystalS_f_absorb, EcalMaterial, "ecalCrystalL_f_absorb", 0, 0, 0);



  std::cout<<"making active with width "<<fiber_diameter<<std::endl;
  G4Box *ecalCrystalS_f_fiber = new G4Box("ecalCrystalS_f_fiber", 0.5*ecal_front_face, 0.5 * ecal_rear_face, 0.5 *fiber_diameter);

  G4LogicalVolume *ecalCrystalL_f_fiber_scinti = new G4LogicalVolume(ecalCrystalS_f_fiber, ScintiMaterial, "ecalCrystalL_f_fiber_scinti", 0, 0, 0);
  G4LogicalVolume *ecalCrystalL_f_fiber_cherenc = new G4LogicalVolume(ecalCrystalS_f_fiber, CherencMaterial, "ecalCrystalL_f_fiber_cherenc", 0, 0, 0);



  //


  // ECAL physical placement
  G4VPhysicalVolume *ecalCrystalP_f[NECAL_CRYST];
  G4VPhysicalVolume *ecalCrystalP_f_fiber_scinti[NECAL_CRYST];
  G4VPhysicalVolume *ecalCrystalP_f_fiber_cherenc[NECAL_CRYST];


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

      z_pos_f[iCrystal] = iCrystal*( hole_diameter+ fiber_diameter*2);

      std::cout<<"iCrystal "<<iCrystal<<" iX iY "<<iX<<" "<<iY<<" z "<<z_pos_f[iCrystal]<<std::endl;


      sprintf(name, "ecalCrystalP_f_absorb_%d", iCrystal);
      ecalCrystalP_f[iCrystal] = new G4PVPlacement(0, G4ThreeVector(0, 0, z_pos_f[iCrystal]+0.5*hole_diameter), ecalCrystalL_f_absorb, name, worldLV, false, 0);

      sprintf(name, "ecalCrystalP_f_fiber_scinti_%d", iCrystal);
      ecalCrystalP_f_fiber_scinti[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(0, 0,  z_pos_f[iCrystal]+hole_diameter+0.5*fiber_diameter), ecalCrystalL_f_fiber_scinti, name, worldLV, false, 0);

      
      sprintf(name, "ecalCrystalP_f_fiber_cherenc_%d", iCrystal);
      ecalCrystalP_f_fiber_cherenc[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(0, 0, z_pos_f[iCrystal]+hole_diameter+1.5*fiber_diameter), ecalCrystalL_f_fiber_cherenc, name, worldLV, false, 0);


      //




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







  G4VisAttributes *VisFiber = new G4VisAttributes(red);
  VisFiber->SetVisibility(true);
  VisFiber->SetForceWireframe(true);
  ecalCrystalL_f_fiber_scinti->SetVisAttributes(VisFiber);
  ecalCrystalL_f_fiber_cherenc->SetVisAttributes(VisFiber);



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


  /************************************************************************************/


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
