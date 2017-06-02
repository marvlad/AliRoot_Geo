// Includes for compilation of this macro
#if !defined (__CINT__) || defined (__MAKECINT__)
  #include "TMCProcess.h"
  #include "TClassTable.h"
  #include "TROOT.h"
  #include "TSystem.h"
  #include "TGeoManager.h"
  #include "TCanvas.h"
  #include "TPaveStats.h"
  #include "TArrow.h"
  // For drawing of tracks
  // #include "TGLClipSet.h"
  #include "TGLClip.h"
  #include "TGLAnnotation.h"
  #include "TGLSAViewer.h"
  #include "TGFrame.h"  // for TGMainFrame.h
  #include "TPolyLine3D.h"
  #include "TPolyMarker3D.h"
  #include "TView.h"
  #include "TView3D.h"
  #include "TGLViewer.h"
  #include "TVirtualGeoTrack.h"
  #include "TGeoTrack.h"
  #include "TObjArray.h"
  #include "TMacro.h"
  #include "TNamed.h"
  // For geometry
  // defined in TGeoVolume: #include "TGeoVolumeAssembly.h"
  #include "TGeoVolume.h"
  #include "TGeoMedium.h"
  #include "TGeoCompositeShape.h"
  #include "TGeoVolume.h"
  #include "TGeoXtru.h"
  #include "TGeoShape.h"
  #include "TGeoBBox.h"
  #include "TGeoPara.h"
  #include "TGeoTrd1.h"
  #include "TGeoTube.h"
  #include "TGeoCone.h"
  #include "TGeoArb8.h"
  #include "TGeoPcon.h"
  #include "TGeoTorus.h"
  #include "TGeoEltu.h"
  //    
  /*#include "AliRunLoader.h"
  #include "AliStack.h"
  #include "AliFMD.h"
  #include "AliFMDHit.h"
  //#include "FMD/FMDsim/AliFMD.h"
  //#include "FMD/FMDsim/AliFMDHit.h"
  #include "AliVZERO.h" 
  #include "AliVZEROLoader.h" 
  #include "AliVZEROhit.h"
  #include "AliITSgeom.h"
  #include "AliITSInitGeometry.h"
  #include "AliITS.h" 
  #include "AliITShit.h"
  #include "AliITSLoader.h" 
  #include "AliITSmodule.h"
  #include "AliZDCHit.h"
  #include "AliZDC.h"
  #include "AliLoader.h"
  #include "AliADhit.h"
  #include "AliAD.h"
  #include "AliADLoader.h"*/
  #include "TH2I.h"
  #include "TH2F.h"
  #include "TH1F.h"
  #include "TTree.h"
  #include "TClonesArray.h"
  #include "TDatabasePDG.h"
  #include "TParticle.h"
  #include "TFile.h"
  #include "TKey.h"
  #include "TParticlePDG.h"
  /*#include "AliHeader.h"
  #include "STEER/STEER/AliRun.h"
  #include "STEER/STEERBase/AliGenEventHeader.h"
  #include "STEER/STEERBase/AliGenPythiaEventHeader.h"
  #include "STEER/STEERBase/AliGenDPMjetEventHeader.h"
  #include "STEER/STEERBase/AliGenHijingEventHeader.h"*/
  #include "TRandom.h"
  #include "TMath.h"
  #include "TArrayF.h"
 /* #include "AliTrackReference.h"*/
  //
#elif defined (__CINT__)
  #include <TInterpreter.h>
#endif

TGeoRotation * Rx180, * Rz180, * Ry180, * Rx90, * Rx90m, * Ry90m, *Ry90, * Rz90;
enum {kCoordX=0, kCoordY=1, kCoordZ=2};
// void PrintNodeAt(Double_t x, Double_t y, Double_t z);
// void FindNodes(const Double_t x, const Double_t y, const Double_t z, const Double_t tf, const Int_t Coord=2) ; 
// TGeoVolumeAssembly * CreateWarmModule_VMAOI();
void CreateCurvedBundles();
TGeoVolume * CreateWarmModule();
TGeoVolumeAssembly * CreateBLM();
TGeoVolumeAssembly * CreatePipeOvalyzed(const TGeoMedium * mat);
TGeoVolumeAssembly * CreateADAShielding();
TGeoVolumeAssembly * CreateVacuumChamberSupport();
TGeoVolumeAssembly * CreatePmtBoxC() ;
TGeoVolumeAssembly * CreateAD();
TGeoVolume * MakeVolIBeam(const char * volname, const TGeoMedium * mat, const Double_t x, const Double_t y, const Double_t dx, const Double_t dy, const Double_t dz);
Bool_t GetADATwoInstalled() { return kTRUE;}
Bool_t GetADCTwoInstalled() { return kTRUE;}
//_________________________________________________________
TGeoVolume * Make_UProfile(const char * volname, const Double_t L, 
                           const TGeoMedium * medium, const Double_t W, const Double_t H, 
                           const Double_t dw, const Double_t dh)
{
  const Int_t nvertices = 8;
  const Double_t W2 = W*0.5;

  Double_t X[nvertices] = {-W2 , -W2 , -W2 + dw , -W2 +dw , W2 -dw , W2 - dw , W2 , W2 }; 
  Double_t Y[nvertices] = { 0  , H   , H        , dh      , dh     , H       , H  , 0  }; 

  TGeoXtru * shAD_UProfileV = new TGeoXtru(2);
  // shADsuppIBeam->SetNameTitle("shAD_UProfile","shAD_UProfile");
  shAD_UProfileV->DefinePolygon(nvertices, X, Y);
  shAD_UProfileV->DefineSection(0, -0.5*L, 0., 0., 1.0); // index, Z position, offset (x,y) and scale for first section
  shAD_UProfileV->DefineSection(1,  0.5*L, 0., 0., 1.0); // idem, second section

  TGeoVolume * vol = new TGeoVolume(volname, shAD_UProfileV, medium);
  vol->SetLineColor(kAzure-3);
  return vol;
}
//_________________________________________________________
TGeoVolume * Make_UProfileH(const char * volname, TGeoMedium * medium)
{
  const Float_t dw = 0.85;
  const Float_t dh = 0.60;
  TGeoVolumeAssembly * voprofile = new TGeoVolumeAssembly(volname);
  TGeoVolume * volH = Make_UProfile(Form("%s_H", volname), 90.-2*dw-0.01, medium, 10, 5, dw, dh);
  TGeoVolume * volL = gGeoManager->MakeBox(Form("%s_L", volname), medium, 5, 2.5, dw/2.);
  volH->SetLineColor(kAzure-3);
  volL->SetLineColor(kAzure-3);
  voprofile->AddNode(volH, 1);
  voprofile->AddNode(volL, 1, new TGeoTranslation(0, 2.5, -90./2. + dw/2.));
  voprofile->AddNode(volL, 2, new TGeoTranslation(0, 2.5, +90./2. - dw/2.));
  return voprofile;
}
//_________________________________________________________
//
enum ADCPosition_t { kADCInTunnel, kADCInCavern, kADCInBoth};
ADCPosition_t fADCPosition = kADCInTunnel;
Bool_t fADCstruct = kTRUE;
Bool_t fADAstruct = kTRUE;
//
TGeoMaterial *matVacuum;
TGeoMaterial *matAl;
// define some media
TGeoMedium * kMedVacuum;
TGeoMedium * kMedAlu;
TGeoMedium * kMedSteelSh;
TGeoMedium * kMedCopper; // = gGeoManager->GetMedium("AD_Cu_C0");  // Stainless Steel 
TGeoMedium * kMedPVC   ; // = gGeoManager->GetMedium("AD_PVC"  );  // Stainless Steel
TGeoVolume * alice = 0; 

Float_t as[11];
Float_t zs[11];
Float_t fs[11];
Double_t radLength, absLength, a_ad, z_ad;
Int_t    id;

// Taken from ./MWE_EvDisplayAD.C
TGeoManager   * fCustomGeo = 0; //!
TGeoManager   * fFileGeo = 0; //!

class CurvedBundle 
{
  // TODO: 
  // 1) Curved bundle should know:
  //    1.a) about its total length. [DONE]
  //    1.b) If a given a (x,y,z) belong the her.
  //    1.C) The fiber length from a point inside her to both ends of the fiber [DONE]

  private:
    TString fVolName;
    TGeoVolumeAssembly * fFiberBundle; 
    // Rotations and translations:
    TGeoTranslation * fGlobalTranslation; 
    TGeoRotation * fRotationYphi   ; //! 
    TGeoRotation * fRotationX90m   ; //! 
    TGeoRotation * fRotationZtheta ; //! 
    Bool_t fReflected; //!
    // Medium and materials
    TGeoMaterial * fMatVacuum; //!
    TGeoMedium   *    fMedium; //!
    Double_t fBundleRadius;    //! radius of the fiber bundle in cm!!
    Double_t fX1; Double_t fY1; Double_t fZ1; //! First point coordinates
    Double_t fX2; Double_t fY2; Double_t fZ2; //! Second point coordinates
    Double_t fPhiAngle;                       //! Phi angle in degrees
    Double_t fThetaAngle;                     //! Angle between plane of Arc0 and horizontal plane (XZ)
    Double_t fLocalW;                         //! 
    Double_t fMaxAllowedX;                    //! Maximun x in local coordinate system:
    Double_t fLocalC0X, fLocalC0Y, fLocalArc0StartAngle, fLocalArc0EndAngle, fLocalArc0Radius; //! Center of first arc in local plane
    Double_t fLocalC1X, fLocalC1Y, fLocalArc1StartAngle, fLocalArc1EndAngle, fLocalArc1Radius; //! Center of first arc in local plane
    Double_t fLocalC2X, fLocalC2Y, fLocalArc2StartAngle, fLocalArc2EndAngle, fLocalArc2Radius; //! Center of first arc in local plane
    Double_t fLocalC3X, fLocalC3Y, fLocalArc3StartAngle, fLocalArc3EndAngle, fLocalArc3Radius; //! Center of first arc in local plane
    Double_t fStartPointRadius; //!
    Double_t fArc1Beta; //!
    Double_t fLengthTotal;  //!
    Double_t fLengthArc[4]; //! 
    // To compute the early-time
    Double_t fLowerLimitMasterZ_Arc[4]; //!
    Double_t fUpperLimitMasterZ_Arc[4]; //!
    TGeoHMatrix *   fLocalMatrixArc[4]; //!
  public:
  CurvedBundle(const char * name, Double_t x1, Double_t y1, Double_t z1, Double_t x2, Double_t y2, Double_t z2, Double_t r1, Double_t r2, Double_t max_allowed_x, Bool_t reflex) :
    fVolName(name),
    fGlobalTranslation(0),
    fRotationYphi(0),
    fRotationX90m(0),
    fRotationZtheta(0),
    fReflected(reflex),
    fBundleRadius((1.6/2.)/TMath::Sqrt(2.)), //! radius (r) of the half fiber bundle (in cm). 2(r^2)=R^2, With R the radius of the full fiber bundle !!
    fLengthTotal(-1),
    fFiberBundle(new TGeoVolumeAssembly(fVolName))
  {
    SetFirstPoint  ( x1 , y1 , z1 - r1 ); 
    SetLastPoint   ( x2 , y2 , z2 ); 
    SetFirstPointRadius(r1);
    SetLastPointRadius (r2);
    FindPhiAngle();
    SetMaxAllowedX(max_allowed_x);
    FindLocalPlaneW();
    // First compute parameters in 3->2->1 order !
    ComputeArc3Parameters();
    ComputeArc2Parameters();
    ComputeArc1Parameters();
    ComputeArc0Parameters();
    // You can optionaly make a reflection on zx plane:
    if (fReflected) MakeReflectionXZ();
    // After parameters are computed, we create the matrices !
    CreateMatrices();
    CreateMediums();
    // Finally, we create the volumes
    CreateArc3();
    CreateArc2();
    CreateArc1();
    CreateArc0();
    ComputeTotalLength();
    // fCombiTrans  = new TGeoCombiTrans(fX1, fY1, fZ1, fRotation123);
    Print();
  }
  const char * GetName() { return fVolName; } 
  TGeoVolumeAssembly * GetVolume() { return fFiberBundle; } 
  // TGeoCombiTrans * GetCombiTrans() { return fCombiTrans; }
  void SetFirstPoint(Double_t x, Double_t y, Double_t z) { fX1 = x, fY1 = y, fZ1 = z; }
  void SetLastPoint(Double_t x, Double_t y, Double_t z) { fX2 = x, fY2 = y, fZ2 = z; }
  void SetFirstPointRadius(Double_t r) { fStartPointRadius=r;      } 
  void SetLastPointRadius(Double_t r)  { fLocalArc3Radius=r;      } 
  Double_t GetFirstPointRadius()       { return fStartPointRadius; } 
  Double_t GetLastPointRadius()        { return fLocalArc3Radius; } 
  void SetMaxAllowedX(Double_t x) {
    if (x>fX2 || x<fX1) { 
      printf("max_allowed_x (%f) out of valid range [%f;%f]\n", x, fX1, fX2);
      exit(1);
    }

    Double_t cos = TMath::Cos(DegToRad(fPhiAngle));
    if (cos==0) { printf("%s:%d: Error, angle is 90 degrees\n", __FUNCTION__, __LINE__); exit(1); }
    fMaxAllowedX = (x-fX1) / cos;
  }
  void FindPhiAngle() { /* result = atan2 (y,x) * 180 / PI; */ fPhiAngle = - RadToDeg(TMath::ATan2(fZ2-fZ1, fX2-fX1)); }
  void FindLocalPlaneW() { fLocalW = TMath::Sqrt((fX2-fX1)*(fX2-fX1) + (fZ2-fZ1)*(fZ2-fZ1)); }
  Double_t RadToDeg(Double_t rad) { return rad * (180./TMath::Pi()); }
  Double_t DegToRad(Double_t deg) { return deg * (TMath::Pi()/180.); }
  //
  void NegateSwapAngles(Double_t & angle1, Double_t & angle2)
  {
    Double_t aux = angle1;
    angle1 = -angle2;
    angle2 = -aux;
  }
  void MakeReflectionXZ(){
    // Make Y-Centers negative
    fLocalC3Y = -fLocalC3Y;
    fLocalC2Y = -fLocalC2Y;
    fLocalC1Y = -fLocalC1Y;
    //
    NegateSwapAngles(fLocalArc3StartAngle, fLocalArc3EndAngle);
    NegateSwapAngles(fLocalArc2StartAngle, fLocalArc2EndAngle);
    NegateSwapAngles(fLocalArc1StartAngle, fLocalArc1EndAngle);
    // 
    fThetaAngle = - fThetaAngle;
    fY1 = -fY1;

  }
  void CreateMatrices() 
  {
    // TObjArray *   GetListOfMatrices () const
    fRotationX90m      = new TGeoRotation("Rx90m",   0., -90.,   0.);
    fRotationYphi      = new TGeoRotation(Form("fRotationYphi_%s"  , fVolName. Data()), 90, fPhiAngle, -90); 
    fRotationZtheta    = new TGeoRotation(Form("fRotationZtheta_%s", fVolName. Data()), fThetaAngle, 0, 0 ); 
    fGlobalTranslation = new TGeoTranslation(fX1, fY1, fZ1); 
  }
  void CreateMediums()
  {
    // fMatVacuum = new TGeoMaterial("Vacuum", 0,0,0);
    // fMedium = new TGeoMedium("Vacuum10",10, fMatVacuum);
    fMedium = gGeoManager->GetMedium("AD_BC404"); // AD Scin.
  }
  void ComputeArc3Parameters()
  {
    fLocalC3X = fLocalW - GetLastPointRadius();
    fLocalC3Y = fY2-fY1;
    fLocalArc3StartAngle =  0.;
    fLocalArc3EndAngle   =  90.;
    ///////////////////////////////
    // Ry90m = new TGeoRotation("Ry90m",  90., -90., -90.) ;
    // Ry90  = new TGeoRotation("Ry90" ,  90.,  90., -90.) ;

  }
  void ComputeArc2Parameters()
  {
    Double_t r = GetLastPointRadius(); 
    Double_t l = fLocalW - r - fMaxAllowedX;
    fLocalArc2Radius = (l*l + r*r) / (2*r);
    Double_t alpha = RadToDeg(TMath::ATan2(r,l)); // half arc2!

    fLocalC2X = fLocalW - r;
    fLocalC2Y = fY2-fY1 + r - fLocalArc2Radius;
    fLocalArc2StartAngle = 90.; 
    fLocalArc2EndAngle   = 90. + 2.*alpha; 

  }
  void ComputeArc1Parameters()
  {
    // local variables for angles named by greek letters are in radians
    // private class level members for variables are in degrees !

    // Some alias in the local plane
    Double_t alpha = DegToRad(TMath::Abs(180. - fLocalArc2EndAngle));
    Double_t b     = fStartPointRadius;
    Double_t h     = fY2 - fY1;
    Double_t w     = fMaxAllowedX;

    Double_t a     = (0.5 * ( w*w + h*h - b*b )) / ( b + w * TMath::Sin(alpha) + h * TMath::Cos(alpha) );

    Double_t beta  = TMath::ASin(( h - a * TMath::Cos(alpha) ) / ( a + b ));
    Double_t delta = 0.5 * (TMath::PiOver2() - (alpha + beta));
    Double_t R     = a / TMath::Tan(delta) ;
    
    // Center coordinates:
    // fLocalC1X = fX1 - R * TMath::Sin(beta) + b * TMath::Cos(beta);
    // fLocalC1Y = fY1 + R * TMath::Cos(beta) + b * TMath::Sin(beta); 
    fLocalC1X = w - (R * TMath::Cos(alpha));
    fLocalC1Y = h + (R * TMath::Sin(alpha));

    // Arc1 angles:
    fLocalArc1EndAngle   = 360. - RadToDeg(alpha);
    fLocalArc1StartAngle = 360. - RadToDeg(alpha + 2*delta);
    fLocalArc1Radius     = R;

    // Needed for ComputeArc0Parameters:
    fArc1Beta = RadToDeg(beta);

    // printf("****************************************************************\n");
    // printf("****************************************************************\n");
    // printf("//////////////////////// %s ////////////////////////\n", __FUNCTION__);
    // printf("      h: %f |    a: %f |     b: %f | w: %f\n", h, a, b, w);
    // printf("  alpha: %f | beta: %f | delta: %f\n", RadToDeg(alpha), RadToDeg(beta), RadToDeg(delta));
    // printf("  h-acos(alpha): %f\n", h - a * TMath::Cos(alpha));
    // printf("  w-asin(alpha): %f\n", w - a * TMath::Sin(alpha));
    // printf("          a + b: %f\n", a + b);
    // printf("              R: %f\n", R);
    // printf("      fLocalC1X: %f\n", fLocalC1X);
    // printf("      fLocalC1Y: %f\n", fLocalC1Y);
    // printf("****************************************************************\n");
    // printf("****************************************************************\n");

  }
  void ComputeArc0Parameters()
  {
    // local variables for angles named by greek letters are in radians
    // private class level members for variables are in degrees !

    // I need beta from ComputeArc1Parameters
    Double_t beta  = DegToRad(fArc1Beta);
    Double_t phi   = TMath::Abs(DegToRad(fPhiAngle));
    Double_t b     = TMath::Abs(fStartPointRadius);  // This is not a curvature radius!

    Double_t h2    = ( TMath::Power(b,2) * TMath::Power(TMath::Sin(beta),2) ) + TMath::Power((b+b*TMath::Cos(beta)*TMath::Sin(phi)),2);
    Double_t c2    = TMath::Power(b*TMath::Cos(beta)*TMath::Cos(phi),2);
    Double_t  l    = TMath::Sqrt(h2+c2);

    Double_t sigma = TMath::ACos(l/(2*b));
    Double_t theta = TMath::ASin(TMath::Sin(beta)/TMath::Sin(2*sigma));
    
    Double_t R     = b / TMath::Tan(sigma);

    // Coordinates of Arc0 center:
    fLocalC0X =   R;
    fLocalC0Y = - b;
    fLocalArc0StartAngle = 180. - 2 * RadToDeg(sigma);
    fLocalArc0EndAngle   = 180.;
    fLocalArc0Radius     = R;
    fThetaAngle = RadToDeg(theta);

    // printf("****************************************************************\n");
    // printf("****************************************************************\n");
    // printf("            phi: %f\n", RadToDeg(phi));
    // printf("          sigma: %f\n", RadToDeg(sigma));
    // printf("          theta: %f\n", RadToDeg(theta));
    // printf("              l: %f\n", l);
    // printf("              b: %f\n", b);
    // printf("              R: %f\n", R);
    // printf("      fLocalC0X: %f\n", fLocalC0X);
    // printf("      fLocalC0Y: %f\n", fLocalC0Y);
    // printf("****************************************************************\n");
    // printf("****************************************************************\n");
  }
  void CreateArc0()
  {
    // TGeoTorus(Double_t R,Double_t Rmin,Double_t Rmax,Double_t Phi1,Double_t Dphi); 
    printf("============ A ==========\n");
    TString volname = "CurvedBundle_Arc0";
    TGeoTorus    * fArc0Shape = new TGeoTorus(fLocalArc0Radius, 0., fBundleRadius, fLocalArc0StartAngle, fLocalArc0EndAngle - fLocalArc0StartAngle);
    TGeoVolume   *   fArc0Vol = new TGeoVolume(volname, fArc0Shape, fMedium);
    fArc0Vol->SetLineColor(kPink-8);
    printf("============ B ==========\n");
    if (!fFiberBundle) exit(1);
    fArc0Vol -> Print();
    printf("============ C ==========\n");
    // fFiberBundle -> AddNode( fArc0Vol, 1, new TGeoTranslation(fLocalC0X, fLocalC0Y, 0));
    //
    // Combined transformations: Represent a rotation followed by a translation. 
    //
    // For Arc0 the order is different:
    // 1) Translation in the local plane: E.g. (fLocalC3X,fLocalC3Y)
    // 2) Rotation of the local plane (fPhiAngle) or (Rx90 * RzTheta)
    // 3) Translation in global coordinates to P1 (x,y,z) position
    //
    TGeoTranslation * fLocalArc0Tr = new TGeoTranslation(fLocalC0X, fLocalC0Y, 0);
    // TGeoTranslation * fLocalArc0Tr = new TGeoTranslation(0, 0, 0);
    /// -   after=TRUE (default): THIS*ROT
    /// -   after=FALSE         : ROT*THIS
    fRotationZtheta->MultiplyBy(fRotationX90m, 1);
    TGeoCombiTrans  * fLocalArc0Cb = new TGeoCombiTrans (*fGlobalTranslation, *fRotationZtheta);
    // TGeoCombiTrans  * fLocalArc0Cb = new TGeoCombiTrans (*fGlobalTranslation, *fRotationX90m);
    TGeoHMatrix     * fLocalMatrix = new TGeoHMatrix((*fLocalArc0Cb)*(*fLocalArc0Tr));
    //
    fFiberBundle -> AddNode( fArc0Vol, 1, fLocalMatrix);
    // Compute global z-limits for the arc 0
    Double_t localStart[3] = {0.};
    Double_t localEnd  [3] = {0.};
    localStart[0] = fLocalArc0Radius * TMath::Cos( DegToRad(fLocalArc0StartAngle));
    localStart[1] = fLocalArc0Radius * TMath::Sin( DegToRad(fLocalArc0StartAngle));
    localEnd  [0] = fLocalArc0Radius * TMath::Cos( DegToRad(fLocalArc0EndAngle));
    localEnd  [1] = fLocalArc0Radius * TMath::Sin( DegToRad(fLocalArc0EndAngle));
    // printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    // printf("%s:%d Angles: Start= %f End= %f\n", __FUNCTION__, __LINE__, fLocalArc0StartAngle, fLocalArc0EndAngle);
    // printf("%s:%d localStart[]=(%f,%f,%f)\n", __FUNCTION__, __LINE__, localStart[0], localStart[1], localStart[2]);
    // printf("%s:%d localEnd  []=(%f,%f,%f)\n", __FUNCTION__, __LINE__, localEnd  [0], localEnd  [1], localEnd  [2]);

    Double_t masterStart[3] = {0.};
    Double_t masterEnd  [3] = {0.};
    fLocalMatrix->LocalToMaster(localStart, masterStart);
    fLocalMatrix->LocalToMaster(localEnd  , masterEnd  );
    // printf("%s:%d masterStart[]=(%f,%f,%f)\n", __FUNCTION__, __LINE__, masterStart[0], masterStart[1], masterStart[2]);
    // printf("%s:%d masterEnd  []=(%f,%f,%f)\n", __FUNCTION__, __LINE__, masterEnd  [0], masterEnd  [1], masterEnd  [2]);
    // Array indexes below should match Arc number
    const Int_t arc=0;
    fLowerLimitMasterZ_Arc[arc]= TMath::Min(masterStart[2], masterEnd[2]);
    fUpperLimitMasterZ_Arc[arc]= TMath::Max(masterStart[2], masterEnd[2]);
    fLocalMatrixArc[arc]=fLocalMatrix;
    // printf("%s:%d fLowerLimitMasterZ_Arc[%d]=%f\n", __FUNCTION__, __LINE__, arc, fLowerLimitMasterZ_Arc[arc]);
    // printf("%s:%d fUpperLimitMasterZ_Arc[%d]=%f\n", __FUNCTION__, __LINE__, arc, fUpperLimitMasterZ_Arc[arc]);
    // printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
  }
  void CreateArc3()
  {
    // TGeoTorus(Double_t R,Double_t Rmin,Double_t Rmax,Double_t Phi1,Double_t Dphi); 
    TString volname = "CurvedBundle_Arc3";
    TGeoTorus    * fArc3Shape = new TGeoTorus(GetLastPointRadius(), 0., fBundleRadius, fLocalArc3StartAngle, fLocalArc3EndAngle - fLocalArc3StartAngle);
    TGeoVolume   *   fArc3Vol = new TGeoVolume(volname, fArc3Shape, fMedium);
    fArc3Vol->SetLineColor(kCyan);
    if (!fFiberBundle) exit(1);
    fArc3Vol     -> Print();
    // fFiberBundle -> AddNode( fArc3Vol, 1, new TGeoTranslation(fLocalC3X, fLocalC3Y, 0));
    //
    // Combined transformations: Represent a rotation followed by a translation. 
    //
    // This is the roder:
    // 1) Translation in the local plane: E.g. (fLocalC3X,fLocalC3Y)
    // 2) Rotation of the local plane (fPhiAngle) or (Rx90 * RzTheta)
    // 3) Translation in global coordinates to P1 (x,y,z) position
    //
    TGeoTranslation * fLocalArc3Tr = new TGeoTranslation(fLocalC3X, fLocalC3Y, 0);
    TGeoCombiTrans  * fLocalArc3Cb = new TGeoCombiTrans (*fGlobalTranslation, *fRotationYphi);
    TGeoHMatrix     * fLocalMatrix = new TGeoHMatrix((*fLocalArc3Cb)*(*fLocalArc3Tr));
    //
    fFiberBundle -> AddNode( fArc3Vol, 1, fLocalMatrix);
    // Compute global z-limits for the arc 3
    Double_t localStart[3] = {0.};
    Double_t localEnd  [3] = {0.};
    localStart[0] = GetLastPointRadius() * TMath::Cos( DegToRad(fLocalArc3StartAngle));
    localStart[1] = GetLastPointRadius() * TMath::Sin( DegToRad(fLocalArc3StartAngle));
    localEnd  [0] = GetLastPointRadius() * TMath::Cos( DegToRad(fLocalArc3EndAngle));
    localEnd  [1] = GetLastPointRadius() * TMath::Sin( DegToRad(fLocalArc3EndAngle));
    // printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    // printf("%s:%d Angles: Start= %f End= %f\n", __FUNCTION__, __LINE__, fLocalArc3StartAngle, fLocalArc3EndAngle);
    // printf("%s:%d localStart[]=(%f,%f,%f)\n", __FUNCTION__, __LINE__, localStart[0], localStart[1], localStart[2]);
    // printf("%s:%d localEnd  []=(%f,%f,%f)\n", __FUNCTION__, __LINE__, localEnd  [0], localEnd  [1], localEnd  [2]);

    Double_t masterStart[3] = {0.};
    Double_t masterEnd  [3] = {0.};
    fLocalMatrix->LocalToMaster(localStart, masterStart);
    fLocalMatrix->LocalToMaster(localEnd  , masterEnd  );
    // printf("%s:%d masterStart[]=(%f,%f,%f)\n", __FUNCTION__, __LINE__, masterStart[0], masterStart[1], masterStart[2]);
    // printf("%s:%d masterEnd  []=(%f,%f,%f)\n", __FUNCTION__, __LINE__, masterEnd  [0], masterEnd  [1], masterEnd  [2]);
    // Array indexes below should match Arc number
    const Int_t arc=3;
    fLowerLimitMasterZ_Arc[arc]= TMath::Min(masterStart[2], masterEnd[2]);
    fUpperLimitMasterZ_Arc[arc]= TMath::Max(masterStart[2], masterEnd[2]);
    fLocalMatrixArc[arc]=fLocalMatrix;
    // printf("%s:%d fLowerLimitMasterZ_Arc[%d]=%f\n", __FUNCTION__, __LINE__, arc, fLowerLimitMasterZ_Arc[arc]);
    // printf("%s:%d fUpperLimitMasterZ_Arc[%d]=%f\n", __FUNCTION__, __LINE__, arc, fUpperLimitMasterZ_Arc[arc]);
    // printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
  }
  void CreateArc2()
  {
    // TGeoTorus(Double_t R,Double_t Rmin,Double_t Rmax,Double_t Phi1,Double_t Dphi); 
    TString volname = "CurvedBundle_Arc2";
    TGeoTorus    * fArc2Shape = new TGeoTorus(fLocalArc2Radius, 0., fBundleRadius, fLocalArc2StartAngle, fLocalArc2EndAngle - fLocalArc2StartAngle);
    TGeoVolume   *   fArc2Vol = new TGeoVolume(volname, fArc2Shape, fMedium);
    fArc2Vol->SetLineColor(kViolet+10);
    if (!fFiberBundle) exit(1);
    fArc2Vol     -> Print();
    // fFiberBundle -> AddNode( fArc2Vol, 1, new TGeoTranslation(fLocalC2X, fLocalC2Y, 0));
    //
    // Combined transformations: Represent a rotation followed by a translation. 
    //
    // This is the order:
    // 1) Translation in the local plane: E.g. (fLocalC3X,fLocalC3Y)
    // 2) Rotation of the local plane (fPhiAngle) or (Rx90 * RzTheta)
    // 3) Translation in global coordinates to P1 (x,y,z) position
    //
    TGeoTranslation * fLocalArc2Tr = new TGeoTranslation(fLocalC2X, fLocalC2Y, 0);
    TGeoCombiTrans  * fLocalArc2Cb = new TGeoCombiTrans (*fGlobalTranslation, *fRotationYphi);
    TGeoHMatrix     * fLocalMatrix = new TGeoHMatrix((*fLocalArc2Cb)*(*fLocalArc2Tr));
    //
    fFiberBundle -> AddNode( fArc2Vol, 1, fLocalMatrix);
    // Compute global z-limits for the arc 2
    Double_t localStart[3] = {0.};
    Double_t localEnd  [3] = {0.};
    localStart[0] = fLocalArc2Radius * TMath::Cos( DegToRad(fLocalArc2StartAngle));
    localStart[1] = fLocalArc2Radius * TMath::Sin( DegToRad(fLocalArc2StartAngle));
    localEnd  [0] = fLocalArc2Radius * TMath::Cos( DegToRad(fLocalArc2EndAngle));
    localEnd  [1] = fLocalArc2Radius * TMath::Sin( DegToRad(fLocalArc2EndAngle));
    // printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    // printf("%s:%d Angles: Start= %f End= %f\n", __FUNCTION__, __LINE__, fLocalArc2StartAngle, fLocalArc2EndAngle);
    // printf("%s:%d localStart[]=(%f,%f,%f)\n", __FUNCTION__, __LINE__, localStart[0], localStart[1], localStart[2]);
    // printf("%s:%d localEnd  []=(%f,%f,%f)\n", __FUNCTION__, __LINE__, localEnd  [0], localEnd  [1], localEnd  [2]);

    Double_t masterStart[3] = {0.};
    Double_t masterEnd  [3] = {0.};
    fLocalMatrix->LocalToMaster(localStart, masterStart);
    fLocalMatrix->LocalToMaster(localEnd  , masterEnd  );
    // printf("%s:%d masterStart[]=(%f,%f,%f)\n", __FUNCTION__, __LINE__, masterStart[0], masterStart[1], masterStart[2]);
    // printf("%s:%d masterEnd  []=(%f,%f,%f)\n", __FUNCTION__, __LINE__, masterEnd  [0], masterEnd  [1], masterEnd  [2]);
    // Array indexes below should match Arc number
    const Int_t arc=2;
    fLowerLimitMasterZ_Arc[arc]= TMath::Min(masterStart[2], masterEnd[2]);
    fUpperLimitMasterZ_Arc[arc]= TMath::Max(masterStart[2], masterEnd[2]);
    fLocalMatrixArc[arc]=fLocalMatrix;
    // printf("%s:%d fLowerLimitMasterZ_Arc[%d]=%f\n", __FUNCTION__, __LINE__, arc, fLowerLimitMasterZ_Arc[arc]);
    // printf("%s:%d fUpperLimitMasterZ_Arc[%d]=%f\n", __FUNCTION__, __LINE__, arc, fUpperLimitMasterZ_Arc[arc]);
    // printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
  }
  void CreateArc1()
  {
    // TGeoTorus(Double_t R,Double_t Rmin,Double_t Rmax,Double_t Phi1,Double_t Dphi); 
    TString volname = "CurvedBundle_Arc1";
    TGeoTorus    * fArc1Shape = new TGeoTorus(fLocalArc1Radius, 0., fBundleRadius, fLocalArc1StartAngle, fLocalArc1EndAngle - fLocalArc1StartAngle);
    TGeoVolume   *   fArc1Vol = new TGeoVolume(volname, fArc1Shape, fMedium);
    fArc1Vol->SetLineColor(kOrange);
    if (!fFiberBundle) exit(1);
    fArc1Vol -> Print();
    // fFiberBundle -> AddNode( fArc1Vol, 1, new TGeoTranslation(fLocalC1X, fLocalC1Y, 0));
    //
    // Combined transformations: Represent a rotation followed by a translation. 
    //
    // This is the order:
    // 1) Translation in the local plane: E.g. (fLocalC3X,fLocalC3Y)
    // 2) Rotation of the local plane (fPhiAngle) or (Rx90 * RzTheta)
    // 3) Translation in global coordinates to P1 (x,y,z) position
    //
    TGeoTranslation * fLocalArc1Tr = new TGeoTranslation(fLocalC1X, fLocalC1Y, 0);
    TGeoCombiTrans  * fLocalArc1Cb = new TGeoCombiTrans (*fGlobalTranslation, *fRotationYphi);
    TGeoHMatrix     * fLocalMatrix = new TGeoHMatrix((*fLocalArc1Cb)*(*fLocalArc1Tr));
    //
    fFiberBundle -> AddNode( fArc1Vol, 1, fLocalMatrix);
    // Compute global z-limits for the arc 1
    Double_t localStart[3] = {0.};
    Double_t localEnd  [3] = {0.};
    localStart[0] = fLocalArc1Radius * TMath::Cos( DegToRad(fLocalArc1StartAngle));
    localStart[1] = fLocalArc1Radius * TMath::Sin( DegToRad(fLocalArc1StartAngle));
    localEnd  [0] = fLocalArc1Radius * TMath::Cos( DegToRad(fLocalArc1EndAngle));
    localEnd  [1] = fLocalArc1Radius * TMath::Sin( DegToRad(fLocalArc1EndAngle));
    // printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    // printf("%s:%d Angles: Start= %f End= %f\n", __FUNCTION__, __LINE__, fLocalArc1StartAngle, fLocalArc1EndAngle);
    // printf("%s:%d localStart[]=(%f,%f,%f)\n", __FUNCTION__, __LINE__, localStart[0], localStart[1], localStart[2]);
    // printf("%s:%d localEnd  []=(%f,%f,%f)\n", __FUNCTION__, __LINE__, localEnd  [0], localEnd  [1], localEnd  [2]);

    Double_t masterStart[3] = {0.};
    Double_t masterEnd  [3] = {0.};
    fLocalMatrix->LocalToMaster(localStart, masterStart);
    fLocalMatrix->LocalToMaster(localEnd  , masterEnd  );
    // printf("%s:%d masterStart[]=(%f,%f,%f)\n", __FUNCTION__, __LINE__, masterStart[0], masterStart[1], masterStart[2]);
    // printf("%s:%d masterEnd  []=(%f,%f,%f)\n", __FUNCTION__, __LINE__, masterEnd  [0], masterEnd  [1], masterEnd  [2]);
    // Array indexes below should match Arc number
    const Int_t arc=1;
    fLowerLimitMasterZ_Arc[arc]= TMath::Min(masterStart[2], masterEnd[2]);
    fUpperLimitMasterZ_Arc[arc]= TMath::Max(masterStart[2], masterEnd[2]);
    fLocalMatrixArc[arc]=fLocalMatrix;
    // printf("%s:%d fLowerLimitMasterZ_Arc[%d]=%f\n", __FUNCTION__, __LINE__, arc, fLowerLimitMasterZ_Arc[arc]);
    // printf("%s:%d fUpperLimitMasterZ_Arc[%d]=%f\n", __FUNCTION__, __LINE__, arc, fUpperLimitMasterZ_Arc[arc]);
    // printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
  }

  Double_t ComputeTotalLength()
  {
    fLengthArc[0] = fLocalArc0Radius * DegToRad(TMath::Abs(fLocalArc0EndAngle - fLocalArc0StartAngle));
    fLengthArc[1] = fLocalArc1Radius * DegToRad(TMath::Abs(fLocalArc1EndAngle - fLocalArc1StartAngle));
    fLengthArc[2] = fLocalArc2Radius * DegToRad(TMath::Abs(fLocalArc2EndAngle - fLocalArc2StartAngle));
    fLengthArc[3] = fLocalArc3Radius * DegToRad(TMath::Abs(fLocalArc3EndAngle - fLocalArc3StartAngle));

    return fLengthTotal = fLengthArc[0] + fLengthArc[1] + fLengthArc[2] + fLengthArc[3];
  }
  Double_t GetTotalLength() { return (fLengthTotal>0)?fLengthTotal:ComputeTotalLength();}
  Int_t GetArcIndexFromZ(Double_t & z)
  {
    // static Bool_t init_array = 0;
    // static Double_t array[5] = {0};
    // if (!init_array) {
    //   array[0] = fLowerLimitMasterZ_Arc[0];
    //   array[4] = fLowerLimitMasterZ_Arc[3];
    //   for (Int_t i=0; i<4; i++)
    //   {
    //     Double_t z1=fUpperLimitMasterZ_Arc[i  ];
    //     Double_t z2=fLowerLimitMasterZ_Arc[i+1];
    //     array[i+1]=(z1+z2)/2.;
    //     printf("@@@ z1: %8.2f z2: %8.2f array: %8.2f\n", z1, z2, array[i+1]); 
    //   }
    //   init_array=1;
    // }
    // // remenber we are in the negative z-axis
    // if (z>array[0]) { printf("@@@ %s:%d TOO-LOW  | %f > %f\n", __FUNCTION__, __LINE__, z, array[0]); return -1; }
    // if (z<array[4]) { printf("@@@ %s:%d TOO-HIGH | %f < %f\n", __FUNCTION__, __LINE__, z, array[4]); return -1; }

    // for (Int_t i=0; i<4; i++)
    // {
    //   if ((z>=array[i]) && (z<array[i+1])) return i;
    //   
    // }
    // return -1;
    for (Int_t i=0; i<4; i++)
    {
      if ((z>=fLowerLimitMasterZ_Arc[i]) && (z<fUpperLimitMasterZ_Arc[i])) return i;
    }
    // If still here, print debugging info and try to find the closest arc
    for (Int_t i=0; i<4; i++)
    {
      // printf("@@@ %s:%d CHECKING | [ %+8.2f ? %+8.2f ? %+8.2f ] \n", __FUNCTION__, __LINE__, fLowerLimitMasterZ_Arc[i], z, fUpperLimitMasterZ_Arc[i] );
      if ((z>=fLowerLimitMasterZ_Arc[i]) && (z<fUpperLimitMasterZ_Arc[i])) return i;
    }
    if (TMath::Abs(z-fLowerLimitMasterZ_Arc[3])<fBundleRadius+0.1) { z=fLowerLimitMasterZ_Arc[3]; return 3; }
    if (TMath::Abs(z-fLowerLimitMasterZ_Arc[0])<fBundleRadius+0.1) { z=fLowerLimitMasterZ_Arc[0]; return 0; }
    return -1;
  }
  Double_t GetHitLengthInArc(Int_t arc, Double_t x, Double_t y, Double_t z)
  {
    // Transform (x,y,z) to local Arc plane
    Double_t master[3] = {x,y,z};
    Double_t local [3] = {0.};
    fLocalMatrixArc[arc]->MasterToLocal(master, local);
    Double_t R, startangle, endangle, CX, CY;
    switch (arc)
    {
      case 0: R = fLocalArc0Radius; startangle = fLocalArc0StartAngle; endangle = fLocalArc0EndAngle; CX = fLocalC0X; CY = fLocalC0Y; break;
      case 1: R = fLocalArc1Radius; startangle = fLocalArc1StartAngle; endangle = fLocalArc1EndAngle; CX = fLocalC1X; CY = fLocalC1Y; break;
      case 2: R = fLocalArc2Radius; startangle = fLocalArc2StartAngle; endangle = fLocalArc2EndAngle; CX = fLocalC2X; CY = fLocalC2Y; break;
      case 3: R = fLocalArc3Radius; startangle = fLocalArc3StartAngle; endangle = fLocalArc3EndAngle; CX = fLocalC3X; CY = fLocalC3Y; break;
    }
    Double_t dx  = local[0]; // - CX;
    Double_t dy  = local[1]; // - CY;
    Double_t phi = TMath::ATan2(dy, dx);
    phi=RadToDeg(phi);
    Int_t n = 0;
    while (phi<startangle && n++<10) phi+=180.;
    n=0;
    while (phi>endangle   && n++<10) phi-=180;

    // for Arc1 start angle is closer to IP (lower absolute value of z)
    Double_t dangle = (arc==1)?startangle-phi:endangle-phi;
    dangle=TMath::Abs(dangle);
    
    // printf("###########################################################\n");
    // printf("%s:%d: Arc: %d passed-x: %f passed-y: %f passed-z: %f\n", __FUNCTION__, __LINE__, arc, x, y, z);
    // printf("%s:%d: Arc: %d local-x: %f local-y: %f local-z: %f\n", __FUNCTION__, __LINE__, arc, local[0], local[1], local[2]);
    // printf("%s:%d: Arc: %d Star: %f Hit: %f End: %f\n", __FUNCTION__, __LINE__, arc, startangle, phi, endangle);
    // printf("###########################################################\n");

    dangle=DegToRad(dangle);
    return R*dangle;
  }
  Double_t ComputeLengthToWLSBar(Double_t x, Double_t y, Double_t z)
  { 
    // This finds the length from the WLS bar end to the hit.

    // Find out which arc
    Int_t arc = GetArcIndexFromZ(z);

    Double_t lengthPmt2Hit = 0;
    // Transform global (x,y,z) coordinate to local coordinates of the arc.
    // With x-local and y-local and the center (x-c,y-c) of the arc calculate the angle.
    // Withe the angle, local-radius, local-start and local-end angles of the arc, calculate the length
    lengthPmt2Hit += fLengthArc[arc] - GetHitLengthInArc(arc,x,y,z);
    // Add the lenght of the other arcs
    for (Int_t i=3; i>arc; i--)
    {
      lengthPmt2Hit+=fLengthArc[i];
      // printf("%s:%d length: %f\n", __FUNCTION__, __LINE__, lengthPmt2Hit); //dbg-ecv
    }
    //
    return lengthPmt2Hit;
  }
  Double_t ComputeLength(Double_t x, Double_t y, Double_t z)
  { 
    // This finds the length from the position of hit to the end oposite to the WLS bar.
    // Compare to: ComputeLengthToWLSBar()

    // Find out which arc
    Int_t arc = GetArcIndexFromZ(z);

    Double_t lengthPmt2Hit = 0;
    // Transform global (x,y,z) coordinate to local coordinates of the arc.
    // With x-local and y-local and the center (x-c,y-c) of the arc calculate the angle.
    // Withe the angle, local-radius, local-start and local-end angles of the arc, calculate the length
    lengthPmt2Hit += GetHitLengthInArc(arc,x,y,z);
    // Add the lenght of the other arcs
    for (Int_t i=0; i<arc; i++)
    {
      lengthPmt2Hit+=fLengthArc[i];
    }
    //
    return lengthPmt2Hit;
  }
  void TestArc(Int_t arc)
  {
    Double_t R, startangle, endangle, CX, CY;
    switch (arc)
    {
      case 0: R = fLocalArc0Radius; startangle = fLocalArc0StartAngle; endangle = fLocalArc0EndAngle; CX = fLocalC0X; CY = fLocalC0Y; break;
      case 1: R = fLocalArc1Radius; startangle = fLocalArc1StartAngle; endangle = fLocalArc1EndAngle; CX = fLocalC1X; CY = fLocalC1Y; break;
      case 2: R = fLocalArc2Radius; startangle = fLocalArc2StartAngle; endangle = fLocalArc2EndAngle; CX = fLocalC2X; CY = fLocalC2Y; break;
      case 3: R = fLocalArc3Radius; startangle = fLocalArc3StartAngle; endangle = fLocalArc3EndAngle; CX = fLocalC3X; CY = fLocalC3Y; break;
    }
    Double_t radangle = DegToRad(gRandom->Uniform(startangle, endangle));
    Double_t x  = R * TMath::Cos(radangle);
    Double_t y  = R * TMath::Sin(radangle);
    Double_t local[3] = {x,y,0};
    Double_t master[3] = {0};
    fLocalMatrixArc[arc]->LocalToMaster(local, master);
    
    printf("###########################################################\n");
    printf("%s:%d: Arc: %d Star: %f Hiti(Rnd): %f End: %f\n", __FUNCTION__, __LINE__, arc, startangle, RadToDeg(radangle), endangle);
    printf("%s:%d: Arc: %d mx: %f my: %f mz: %f\n", __FUNCTION__, __LINE__, arc, master[0], master[1], master[2]);
    printf("###########################################################\n");

    // Test!!
    ComputeLength(master[0], master[1], master[2]);
  }
  void Print() {
    printf("****************************************************************\n");
    printf("****************************************************************\n");
    printf("   First Point: %f %f %f | R1: %f\n", fX1, fY1, fZ1, GetFirstPointRadius());
    printf("  Second Point: %f %f %f | R2: %f\n", fX2, fY2, fZ2, GetLastPointRadius());
    printf("  Phi Angle is: %f (degrees)\n", fPhiAngle);
    printf("  fMaxAllowedX (in local frame) is: %f\n", fMaxAllowedX);
    printf("  fLocalW (in local frame) is: %f\n", fLocalW);
    printf("  Arc3 parameters: %f %f %f %f R: %f\n", fLocalC3X, fLocalC3Y, fLocalArc3StartAngle, fLocalArc3EndAngle, fLocalArc3Radius);
    printf("  Arc2 parameters: %f %f %f %f R: %f\n", fLocalC2X, fLocalC2Y, fLocalArc2StartAngle, fLocalArc2EndAngle, fLocalArc2Radius);
    printf("  Arc1 parameters: %f %f %f %f R: %f\n", fLocalC1X, fLocalC1Y, fLocalArc1StartAngle, fLocalArc1EndAngle, fLocalArc1Radius);
    printf("  Arc0 parameters: %f %f %f %f R: %f\n", fLocalC0X, fLocalC0Y, fLocalArc0StartAngle, fLocalArc0EndAngle, fLocalArc0Radius);
    printf("  Total Length: %f\n", fLengthTotal);
    for (Int_t i=0; i<4; i++)
    {
      printf("  Length  Arc%d: %f\n", i, fLengthArc[i]);
    }

    printf("****************************************************************\n");
    printf("****************************************************************\n");

  }
};

void PrintVolumes()
{
  gGeoManager = fFileGeo;
  TGeoVolume *vol;
  TGeoVolume * topvol = gGeoManager->GetVolume("ALIC"); 
  Int_t nd = topvol->GetNdaughters();

  for (Int_t i=0; i<nd; i++) {
    gGeoManager->RestoreMasterVolume();
    TGeoNode * node;
    node = topvol->GetNode(i);
    vol = node->GetVolume();
    TString title = node->GetName();
    printf("node: %s\n", node->GetName());
  }

  gGeoManager = fCustomGeo;
}
void AddVolumeToCustomGeo(const char * path)
{
  TGeoManager * aux = NULL;
  aux = gGeoManager;
  gGeoManager = fFileGeo;
  printf("========== AddVolumeFromManager: %s ==========\n", path);
  if (!gGeoManager) exit(1);
  gGeoManager->cd(path);
  TGeoHMatrix * matrix = (TGeoHMatrix *) gGeoManager->GetCurrentMatrix()->MakeClone();
  TGeoNode    * node   = gGeoManager->GetCurrentNode();
  TGeoVolume  * vol    = gGeoManager->GetCurrentNode()->GetVolume()->CloneVolume();
  gGeoManager = fCustomGeo;
  alice->AddNode(vol  , 1, matrix  );

}

void ImportGeoFromFile() 
{
  TGeoManager * aux = NULL;
  if (!fFileGeo) {
    aux = gGeoManager;
    gGeoManager = 0;
    TGeoManager::Import(Form("./geometry.root"));
    fFileGeo = gGeoManager;
  }
  gGeoManager = aux;
}
// 
// NewAD2016Aug_Marvin_2.C
//
void NewAD2016Aug_Marvin_3()
{
  // gStyle->SetCanvasPreferGL(kTRUE);
  gSystem->Load("libGeom");

  // if (gGeoManager) gGeoManager = 0;
  fCustomGeo = new TGeoManager("NewAD", "Test for New AD geometry");
  //
  // Make the top container volume
  //
  alice = gGeoManager -> MakeBox("TOP", kMedVacuum, 400., 400., 2000.);
  // fCustomTop = alice;
  gGeoManager -> SetTopVolume(alice);

  matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
  matAl = new TGeoMaterial("Al", 26.98,13,2.7);
  kMedVacuum = new TGeoMedium("Vacuum",1, matVacuum);
  kMedAlu = new TGeoMedium("Aluminium",2, matAl);
  kMedSteelSh  = kMedAlu; // TODO: Change this!!
  kMedCopper   = kMedAlu; // TODO: Change this!! // = gGeoManager->GetMedium("AD_Cu_C0");  // Stainless Steel 
  kMedPVC      = kMedAlu; // TODO: Change this!! // = gGeoManager->GetMedium("AD_PVC"  );  // Stainless Steel
  TGeoMedium * medADASci     = gGeoManager->GetMedium("AD_NE102"); // AD Scin.

 
  if (1) {
    Rx90m = new TGeoRotation("Rx90m",   0., -90.,   0.) ;
    Rx90  = new TGeoRotation("Rx90" ,   0.,  90.,   0.) ;
    Rx180 = new TGeoRotation("Rx180",   0., 180.,   0.) ;  //   4    |   1
    Rz180 = new TGeoRotation("Rz180", 180.,   0.,   0.) ;  // --------------->  x  
    Rz90  = new TGeoRotation("Rz90" ,  90.,   0.,   0.) ;
    Ry180 = new TGeoRotation("Ry180", 180., 180.,   0.) ;  //   3    |   2
    Ry90m = new TGeoRotation("Ry90m",  90., -90., -90.) ;
    Ry90  = new TGeoRotation("Ry90" ,  90.,  90., -90.) ;
    TGeoVolume * voTest = CreateWarmModule();
    // TGeoVolume * voTest = CreateOldADA();
    // TGeoVolume * voTest = CreatePump();
    voTest -> SetTransparency(0);
    // voTest->SetLineColor(kGray);
    alice->AddNode(voTest      , 1, new TGeoTranslation(0., 0., 0.));
    // gGeoManager->DefaultColors();
    gGeoManager->CloseGeometry();
    gGeoManager->SetVisLevel(50);
    gGeoManager->SetVisOption(0);
    alice->Draw("ogl");
    return;
  }
  // // CreateCurvedBundles();
  CreateAD();

  ImportGeoFromFile();
  // PrintVolumes();
  // AddVolumeToCustomGeo("/ALIC_1/ZEM_1"      ); 
  // AddVolumeToCustomGeo("/ALIC_1/ZEM_2"      ); 
  // AddVolumeToCustomGeo("/ALIC_1/RB24_1"     ); 
  // AddVolumeToCustomGeo("/ALIC_1/DCM0_1"     ); 
  //
  gGeoManager->CloseGeometry();
  gGeoManager->SetVisLevel(20);
  gGeoManager->SetVisOption(0);
  alice->Draw("ogl");
  return;
  // for (Int_t i=0; i<10; i++)
  // {
  //   printf("=========== test %3d ============\n", i);
  //   fBottomBundles[2]->TestArc(1);
  // }
  
  // for (Int_t i=2; i<fNBundles; i++)
  // {
  //   Double_t L = fTopBundles[i]->GetTotalLength();
  //   printf("fTopBundles[%d]: L: %6.3f | 250-L: %6.3f\n", i, L, 250. - L);
  // }
  // for (Int_t i=2; i<fNBundles; i++)
  // {
  //   Double_t L = fBottomBundles[i]->GetTotalLength();
  //   printf("fBottomBundles[%d]: L: %6.3f | 250-L: %6.3f\n", i, L, 250. - L);
  // }
  //
  // CurvedBundle bundle1("Bundle1",  -95,0,-1920, 20,25,-1955, 20, 20, -60);
  // CurvedBundle bundle2("Bundle2",  -95,0,-1920,  2,25,-1955, 20, 20, -60);
  // CurvedBundle bundle3("Bundle3", -100,0,-1920, -2,25,-1955, 20, 20, -60);
  // CurvedBundle bundle4("Bundle4", -100,0,-1920,-20,25,-1955, 20, 20, -60);
  // CurvedBundle bundle1r("Bundle1r",  -95,0,-1920, 20,25,-1955, 20, 20, -60, kTRUE );
  // alice -> AddNode( bundle1.  GetVolume(), 1);
  // alice -> AddNode( bundle2. GetVolume(), 2);
  // alice -> AddNode( bundle3. GetVolume(), 3);
  // alice -> AddNode( bundle4. GetVolume(), 4);
  // alice -> AddNode( bundle1. GetVolume(), 1, bundle1. GetCombiTrans());



  // Add volumes from geometry file:
  ImportGeoFromFile();
  PrintVolumes();
  gGeoManager->SetVisOption(0);
  AddVolumeToCustomGeo("/ALIC_1/DCM0_1"     ); 
  AddVolumeToCustomGeo("/ALIC_1/YSAA_1"     ); 
  AddVolumeToCustomGeo("/ALIC_1/RB26Pipe_1" ); 
  AddVolumeToCustomGeo("/ALIC_1/RB24_1"     ); 
  AddVolumeToCustomGeo("/ALIC_1/HALL_1"     ); 
  AddVolumeToCustomGeo("/ALIC_1/VZERO_1"     ); 
  // alice -> Draw("ogl");
  // return;
  //
  // Define some materials
  //
  matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
  matAl = new TGeoMaterial("Al", 26.98,13,2.7);

  // **************************************************************** 
  //     Defines tracking media parameters.                         
  //     Les valeurs sont commentees pour laisser le defaut         
  //     a GEANT (version 3-21, page CONS200), f.m.                 
  Float_t epsil, stmin, tmaxfd, deemax, stemax;                     
  epsil  = .001 ; // Tracking precision,
  stemax = -1.  ; // Maximum displacement for multiple scat
  tmaxfd = -20. ; // Maximum angle due to field deflection
  deemax = -.3  ; // Maximum fractional energy loss, DLS
  stmin  = -.8  ; 
  // **************************************************************** 
  
  //
  // Parameters  for AD scintillator: NE-102 (BC400)
  //
  // NE-102, has the following properties : (from internet, need better reference)
  //    Density : ca. 1.03 g/cm3
  //    Electrons/cm3: 3.39 x 10^23
  //    H atoms/cm3: 5.28 x 10^22
  //    C atoms/cm3: 4.78 x 10^22
  //    Ratio of H to C : 1.104 .
  //    wavelength of emission : ~4.23 nm.
  //    Decay time : ~2.4 ns.
  //    Luminescent efficiency : typically 18% of NaI(Tl)
  //    Photons/MeV: 2.5 x 10^4 
  //
  // H                // C 
  as[0] = 1.00794;    as[1] = 12.011;
  zs[0] = 1.;         zs[1] = 6.;
  fs[0] = 5.23;       fs[1] = 4.74;
  Double_t density = 1.032;
  id      = 1;
  //  AliMixture( id, "NE102", as, zs, density, -2, fs );
  //  AliMedium( id, "NE102", id, 1, fieldType, maxField, maxBending, maxStepSize,
  //             maxEnergyLoss, precision, minStepSize );


  //
  // Parameters for lightGuide:  
  //     TODO check material 
  // Should be Poly(methyl methacrylate) (PMMA) acrylic 
  // (C5O2H8)n 
  // Density  1.18 g/cm3
  // Mixture PMMA    Aeff=12.3994 Zeff=6.23653 rho=1.18 radlen=34.0677 intlen=63.3073
  // Element #0 : C  Z=  6.00 A= 12.01 w= 0.600 natoms=5
  // Element #1 : H  Z=  1.00 A=  1.01 w= 0.081 natoms=8
  // Element #2 : O  Z=  8.00 A= 16.00 w= 0.320 natoms=2

  // Carbon          Hydrogen          Oxygen
  as[0] = 12.0107;   as[1] = 1.00794;  as[2] = 15.9994;
  zs[0] = 6.;        zs[1] = 1.;       zs[2] = 8.;
  fs[0] = 0.60;      fs[1] = 0.081;    fs[2] = 0.32;
  density = 1.18;
  id      = 2;
//  AliMixture( id, "PMMA", as, zs, density, 3, fs );
//  AliMedium( id,"PMMA", id, 1, fieldType, maxField, maxBending, maxStepSize,
//             maxEnergyLoss, precision, minStepSize );



  // //
  // // Define some media
  // //
  // kMedVacuum = new TGeoMedium("Vacuum",1, matVacuum);
  // kMedAlu = new TGeoMedium("Aluminium",2, matAl);
  // kMedSteelSh  = kMedAlu;  // TODO: Change this!!
  // TGeoMedium * medADASci     = gGeoManager->GetMedium("AD_NE102"); // AD Scin.

  // CreateAD();
  // // CreateWarmModule_VMAOI();
  // printf("==============DEBUG=============\n");
  // return ;
  gGeoManager->DefaultColors();
  gGeoManager->CloseGeometry();
  gGeoManager->SetVisLevel(10);
  gGeoManager->SetVisOption(0);
  alice->Draw("ogl");

}
//_____________________________________________________________________________
TGeoCompositeShape * MakeShapeADCpadH(const Double_t W, const Double_t H, const Double_t dz) {
  /////////////////////////////////////////////////////////////////////////////
  ///                ADC pad in the cavern (H shapped hole)                 ///
  /////////////////////////////////////////////////////////////////////////////
  // const Double_t W = kADCCellSide; // Width  of Scintillator pad
  // const Double_t H = kADCCellSide; // Height of Scintillator pad
  // Coordinates of ADC pad vertexes
  Double_t pad0_x [] = { 10., 10.,   W,  W  };
  Double_t pad0_y [] = {  0., 11., 11.,  0. };
  Double_t pad1_x [] = { 15., 15.,   W,  W  };
  Double_t pad1_y [] = { 11., 15., 15., 11. };
  Double_t pad2_x [] = {  0.,  0.,   W,  W  };
  Double_t pad2_y [] = { 15.,  H ,   H, 15. };
  TGeoArb8 * shADCpad0H = new TGeoArb8("shADCpad0H", dz/2.);
  TGeoArb8 * shADCpad1H = new TGeoArb8("shADCpad1H", dz/2.);
  TGeoArb8 * shADCpad2H = new TGeoArb8("shADCpad2H", dz/2.);
  for (Int_t i=0; i<4; i++) {
    // -dz
    shADCpad0H -> SetVertex(i,   pad0_x[i], pad0_y[i]);
    shADCpad1H -> SetVertex(i,   pad1_x[i], pad1_y[i]);
    shADCpad2H -> SetVertex(i,   pad2_x[i], pad2_y[i]);
    // +dz
    shADCpad0H -> SetVertex(i+4, pad0_x[i], pad0_y[i]);
    shADCpad1H -> SetVertex(i+4, pad1_x[i], pad1_y[i]);
    shADCpad2H -> SetVertex(i+4, pad2_x[i], pad2_y[i]);
  }
  return new TGeoCompositeShape("shADCpadH","shADCpad0H + shADCpad1H + shADCpad2H");
}
//_____________________________________________________________________________
TGeoCompositeShape * MakeFiberBundle(const char * shName, const Double_t xWLS, const Double_t xPMT, const Double_t yWLS, const Double_t yPMT, const Double_t Z1)
{
  const Double_t fFiberBundleDx =  0.40;
  const Double_t fFiberBundleDz =  2.54;
  const Double_t R = fFiberBundleDx; 

  Double_t fFiberBundleHT = yPMT-yWLS; // Total Height
  Double_t fFiberBundleHW = 5.00; // Height of Fiber Bundle close to WLS
  Double_t fFiberBundleHP = 5.00; // Height of Fiber Bundle close to PMT
  Bool_t cond  = 0;
  Bool_t condY = 0;
  if (xPMT>xWLS)    cond =1;
  else              cond =0;
  if      (yWLS>0 && yPMT>0) condY = 1;
  else if (yWLS<0 && yPMT<0) condY = 0;
  else {
    printf("I'm not prepared to handle this case!!\nQuiting!!\n");
    exit(1);
  }
  // 
  Double_t W = 0;
  Double_t pL, pR;
  if (cond) {
    pL = xWLS + 0.5*R; 
    pR = xPMT -     R; 
  } else {
    pL = xPMT +     R; 
    pR = xWLS - 0.5*R; 
  }
  Double_t Up, Dw;
  Double_t CyWLS, CyPMT;
  if (condY) {
    // Upper side
    Up = yPMT - TMath::Abs(fFiberBundleHP);
    Dw = yWLS + TMath::Abs(fFiberBundleHW);
    // Center of bendings
    CyWLS = yWLS + TMath::Abs(fFiberBundleHW);
    CyPMT = yPMT - TMath::Abs(fFiberBundleHP);
  } else {
    // Down side
    Up = yWLS - TMath::Abs(fFiberBundleHW);
    Dw = yPMT + TMath::Abs(fFiberBundleHP);
    // Center of bendings
    CyWLS = yWLS - TMath::Abs(fFiberBundleHW);
    CyPMT = yPMT + TMath::Abs(fFiberBundleHP);
  }
  
  Double_t A = pR -pL;
  Double_t H = Up - Dw;
  if (!condY) {
    fFiberBundleHW *= -1.;
    fFiberBundleHP *= -1.;
  }
  // angles
  Double_t D = TMath::Sqrt(H*H + A*A);
  Double_t L = TMath::Sqrt(D*D - R*R);
  Double_t beta = TMath::ATan2(A, H);
  Double_t gamm = TMath::ATan2(R, L);
  Double_t alph = beta+gamm;
  // displacement of middle section center:
  Double_t Cx = 0;
  Double_t Cy = 0;
  Cy = (Up + Dw)/2.;
  Cx = pL + A/2.;
  Double_t OriginM[3] = {Cx,Cy,Z1};
  Double_t OriginW[3] = {xWLS,yWLS+fFiberBundleHW/2.,Z1};
  Double_t OriginP1[3] = {xPMT-R/2., yPMT-fFiberBundleHP/2.,Z1};
  Double_t OriginP2[3] = {xPMT+R/2., yPMT-fFiberBundleHP/2.,Z1};
  Double_t * OriginP = 0;
  if (cond ) OriginP = OriginP1;
  else       OriginP = OriginP2;
  Double_t beta_deg = beta*TMath::RadToDeg();
  Double_t gamm_deg = gamm*TMath::RadToDeg();
  Double_t alph_deg = alph*TMath::RadToDeg();
  printf("===========================================================\n");
  printf(" beta = %f (deg)\n", beta_deg);
  printf(" gamm = %f (deg)\n", gamm_deg);
  printf(" alph = %f (deg)\n", alph_deg);
  printf("    L = %f  (cm)\n",  L);
  printf("   Cx = %f  (cm)\n", Cx);
  printf("   Cy = %f  (cm)\n", Cy);
  printf("===========================================================\n");
  TString fNameMid     = Form("sh_%s_FiberBundleMiddle"  , shName ); 
  TString fNameWLS     = Form("sh_%s_FiberBundleNearWLS" , shName ); 
  TString fNamePMT     = Form("sh_%s_FiberBundleNearPMT" , shName ); 
  TString fNameBendWLS = Form("sh_%s_FiberBendWLS"       , shName ); 
  TString fNameBendPMT = Form("sh_%s_FiberBendPMT"       , shName ); 

  TGeoBBox * shFiberBundleMiddle  = new TGeoBBox(fNameMid , R/2.              , L/2.              , fFiberBundleDz/2.); 
  TGeoBBox * shFiberBundleNearWLS = new TGeoBBox(fNameWLS , fFiberBundleDx/2. , TMath::Abs(fFiberBundleHW/2.) , fFiberBundleDz/2., OriginW); 
  TGeoBBox * shFiberBundleNearPMT = new TGeoBBox(fNamePMT , fFiberBundleDx/2. , TMath::Abs(fFiberBundleHP/2.) , fFiberBundleDz/2., OriginP); 

  Double_t fOriginWLS[3] = {0.};
  Double_t fOriginPMT[3] = {0.};
  //
  // WLS Bending Setup
  //
  fOriginWLS[1] = CyWLS;
  fOriginWLS[2] = Z1;
  if (cond ) fOriginWLS[0] = xWLS + R/2.;
  else       fOriginWLS[0] = xWLS - R/2.;
  if (condY) {
    if (cond ) new TGeoTubeSeg(fNameBendWLS, 0, R, fFiberBundleDz/2., 180. - alph_deg,     180.);
    else       new TGeoTubeSeg(fNameBendWLS, 0, R, fFiberBundleDz/2.,              0., alph_deg);
  } else {
    if (cond ) new TGeoTubeSeg(fNameBendWLS, 0, R, fFiberBundleDz/2., 180.,      180.+ alph_deg);
    else       new TGeoTubeSeg(fNameBendWLS, 0, R, fFiberBundleDz/2., 360. - alph_deg,     360.);
  }
  TString fNameBendWLStra = "Tra_" + fNameBendWLS;
  (new TGeoTranslation(fNameBendWLStra, fOriginWLS[0], fOriginWLS[1], fOriginWLS[2]))->RegisterYourself();
  //
  // PMT Bending Setup
  //
  fOriginPMT[1] = CyPMT;
  fOriginPMT[2] = Z1;
  if (cond ) fOriginPMT[0] = xPMT - R;
  else       fOriginPMT[0] = xPMT + R;
  if (condY) {
    if (cond ) new TGeoTubeSeg(fNameBendPMT, 0, R, fFiberBundleDz/2., 360. - alph_deg,           360.);
    else       new TGeoTubeSeg(fNameBendPMT, 0, R, fFiberBundleDz/2.,            180., 180 + alph_deg);
  } else {
    if (cond ) new TGeoTubeSeg(fNameBendPMT, 0, R, fFiberBundleDz/2.,              0.,       alph_deg);
    else       new TGeoTubeSeg(fNameBendPMT, 0, R, fFiberBundleDz/2., 180. - alph_deg,           180.);
  }
  TString fNameBendPMTtra = "Tra_" + fNameBendPMT;
  (new TGeoTranslation(fNameBendPMTtra, fOriginPMT[0], fOriginPMT[1], fOriginPMT[2]))->RegisterYourself();
  //
  // Name for CombiTrans
  TString fNameMidRot = "Rot_" + fNameMid;
  TString fNameMidCom = "Com_" + fNameMid;
  TGeoRotation * fRotM = 0;
  if (condY)
  {
    if (cond ) fRotM = new TGeoRotation(fNameMidRot, -alph_deg, 0., 0.);
    else       fRotM = new TGeoRotation(fNameMidRot,  alph_deg, 0., 0.);
  } else { 
    if (cond ) fRotM = new TGeoRotation(fNameMidRot,  alph_deg, 0., 0.);
    else       fRotM = new TGeoRotation(fNameMidRot, -alph_deg, 0., 0.);
  }
  (new TGeoCombiTrans(fNameMidCom, OriginM[0], OriginM[1], OriginM[2], fRotM))->RegisterYourself();
  //
  TString fExpCompShape = "(";
  fExpCompShape +=         fNameWLS; 
  fExpCompShape += " + " + fNameBendWLS + ":" + fNameBendWLStra ;
  fExpCompShape += " + " + fNameMid     + ":" + fNameMidCom     ;
  fExpCompShape += " + " + fNameBendPMT + ":" + fNameBendPMTtra ;
  fExpCompShape += " + " + fNamePMT     + ")"   ;
  //
  return new TGeoCompositeShape(Form("shFiberBundle_%s", shName), fExpCompShape);

}

//Used later in step manager GetDistanceFromFibersToWLSADA(Int_t ADmodule, Bool_t IsShort, Float_t fHitY)
Double_t fX1FiberShort [8] = {0.};
Double_t fX1FiberLong  [8] = {0.};
Double_t fX2Fiber      [8] = {0.};
Double_t fY1Fiber      [8] = {0.};
Double_t fY2Fiber      [8] = {0.};
// To be moved to AliADv1.h
//
// Private comunication by Arturo Tauro (2014, Apr 23)
// According to last survey measurement done this morning, 
// the C-Side wall is at Z = - 18959mm.
//
Double_t RadToDeg(Double_t rad) { return rad * (180./TMath::Pi()); }
Double_t DegToRad(Double_t deg) { return deg * (TMath::Pi()/180.); }
const Double_t kZwall        = 1895.9 ;  // Aluminium plate z position 
const Double_t kZendAbs      = 1880.75;  // End of CC block absorber
const Double_t kZbegVMAOI    = 1919.2 ;  // Begining of Warm Module
const Double_t kZbegValve    = 1910.7 ;  // Begining of Valve
const Double_t kZbegFrontBar = 1949.1 ;  // Begining of Front Bar
const Double_t kZbegCoil     = 1959.4 ;  // Begining of compensator coil
//_____________________________________________________________________________
TGeoVolumeAssembly * CreatePmtBoxC() 
{
  // Dimensions of the box:
  // Y: 100 cm
  // Z:  70 cm
  // X:  15 cm
  Float_t thickness = 0.5; // cm

  // 
  // Coordinates of the center of the box:
  //
  Float_t xBox = +120. +15./2. +0.01; // The box is touching the shielding block (which goes from x=-120 to 120 cm)
  Float_t yBox = 0;           // At the same level as beam pipe?
  Float_t zBox = kZendAbs - 6.5 - 35.; // The far end of he box is 6.5 cm appart from the far end of shielding block

  TGeoBBox * shBox0 = new TGeoBBox("shAD_PmtBox_C_0", 7.5 + thickness, 50.0 + thickness , 35.0 + thickness );
  TGeoBBox * shBox1 = new TGeoBBox("shAD_PmtBox_C_1", 7.5, 50., 35.);
  Double_t origin[] = {-1,0,2};
  TGeoBBox * shBox2 = new TGeoBBox("shAD_PmtBox_C_Entrance", 3, 49., 35., origin);
  TGeoCompositeShape * shAD_PmtBox_C = new TGeoCompositeShape("shAD_PmtBox_C", "(shAD_PmtBox_C_0-shAD_PmtBox_C_1)-shAD_PmtBox_C_Entrance");

  TGeoVolume * voAD_PmtBox_C = new TGeoVolume("voAD_PmtBox_C", shAD_PmtBox_C, kMedAlu);
  voAD_PmtBox_C->SetLineColor(kGray);
  // 
  // Assemble everything
  //
  TGeoVolumeAssembly * volBox = new TGeoVolumeAssembly("voAD_PMTBOX_C");
  // Make the PM box
  volBox->AddNode(voAD_PmtBox_C, 1, new TGeoTranslation( xBox, yBox, zBox));

  // Make the PMs
  Double_t fPMTdz = 6.;
  TGeoMedium * medADCSci = gGeoManager->GetMedium("AD_BC404");
  TGeoVolume * voADCpmt = new TGeoVolume("voADCpmt", new TGeoTube("shADCpmt", 0., 2.0, fPMTdz/2.), medADCSci);
  voADCpmt -> SetLineColor(kGray);
  Float_t x=xBox;
  // PMT input parameters;
  Float_t a = 20; // Distance (z) from PMT-Box-Wall to top-PMT-center
  Float_t b = 30; // Distance (z) from PMT-Box-Wall to bottom-PMT-center
  Float_t m = 10; // Distance (y) from beam-pipe-center to bottom-PMT-center
  Float_t d =  8; // Distance between PMT-centers
  // PMT calculated parameters:
  Float_t h = TMath::Sqrt((d*3)*(d*3)+(b-a)*(b-a));
  Float_t dzC = (b-a)/3;
  Float_t dyC =     h/3;
  // End of Fiber bundles
  Double_t zBundle   = +1895.;
  Double_t xBundle   = +125.5;
  Double_t yBundle   = 0.;
  Double_t yB[]      = { 7, 5, 3, 1};
  x=xBundle;
  for (Int_t i=0; i<4; i++) 
  {
    Double_t y = m + (dyC * i);
    Double_t z = zBox + 35 - b + (dzC * i);
    Double_t dz =   zBundle - z;
    Double_t yBundle = yB[3-i];
    Double_t dy =   yBundle - y;
    Double_t angle = -(TMath::ATan2(dy,dz));
    Double_t l = TMath::Sqrt(dz*dz+dy*dy)-fPMTdz/2.-1;
    TGeoVolume * voADCfiber = gGeoManager->MakeTube(Form("voAD_FiberPmtC_top_%d", i), medADCSci, 0., 1.6/2., l/2.);
    voADCfiber -> SetLineColor(kCyan);
    printf("@@@@: %s: angle[%d]=%f Length: %f\n", __FUNCTION__, i, RadToDeg(angle), l+1.);
    TGeoRotation * RotXangle = new TGeoRotation(Form("RxADpmtC_top_%d", i),   0., RadToDeg(angle),   0.) ;
    volBox->AddNode(voADCpmt   , i+1 , new TGeoCombiTrans(x , y, z  , RotXangle)); 
    y = yBundle + 0.5*(l+1)*TMath::Sin(angle);
    z = zBundle - 0.5*(l+1)*TMath::Cos(angle);
    volBox->AddNode(voADCfiber , i+1 , new TGeoCombiTrans(x , y, z,  RotXangle)); 
  }
  //
  for (Int_t i=0; i<4; i++) 
  {
    Double_t y = -m - (dyC * i);
    Double_t z = zBox + 35 - b + (dzC * i);
    Double_t dz =   zBundle - z;
    Double_t yBundle = - yB[3-i];
    Double_t dy =   yBundle - y;
    Double_t angle = -(TMath::ATan2(dy,dz));
    Double_t l = TMath::Sqrt(dz*dz+dy*dy)-fPMTdz/2.-1;
    TGeoVolume * voADCfiber = gGeoManager->MakeTube(Form("voAD_FiberPmtC_bot_%d", i), medADCSci, 0., 1.6/2., l/2.);
    voADCfiber -> SetLineColor(kCyan);
    printf("@@@@: %s: angle[%d]=%f Length: %f\n", __FUNCTION__, i, RadToDeg(angle), l+1.);
    TGeoRotation * RotXangle = new TGeoRotation(Form("RxADpmtC_bot_%d", i),   0., RadToDeg(angle),   0.) ;
    volBox->AddNode(voADCpmt   , i+1 , new TGeoCombiTrans(x , y, z  , RotXangle)); 
    y = yBundle + 0.5*(l+1)*TMath::Sin(angle);
    z = zBundle - 0.5*(l+1)*TMath::Cos(angle);
    volBox->AddNode(voADCfiber , i+1 , new TGeoCombiTrans(x , y, z,  RotXangle)); 
  }


  return volBox;

  // fTopBundles[0] = new CurvedBundle     ( "BundleTopRO_l"        , -127 , 7 , fBundleStartZ , 18  , 25 , -1956 , 15 , 22 , -75 , kFALSE ); 
  // fTopBundles[2] = new CurvedBundle     ( "BundleTopRO_s"        , -127 , 5 , fBundleStartZ , 2   , 25 , -1956 , 15 , 18 , -70 , kFALSE ); 
  // fTopBundles[4] = new CurvedBundle     ( "BundleTopLO_s"        , -127 , 3 , fBundleStartZ , -2  , 25 , -1956 , 20 , 15 , -65 , kFALSE ); 
  // fTopBundles[6] = new CurvedBundle     ( "BundleTopLO_l"        , -127 , 1 , fBundleStartZ , -18 , 25 , -1956 , 20 , 10 , -60 , kFALSE ); 

  // fTopBundles[1] = new CurvedBundle     ( "BundleTopRI_l"        , -124 , 7 , fBundleStartZ , 18  , 25 , -1952 , 14 , 22 , -75 , kFALSE ); 
  // fTopBundles[3] = new CurvedBundle     ( "BundleTopRI_s"        , -124 , 5 , fBundleStartZ , 2   , 25 , -1952 , 14 , 18 , -70 , kFALSE ); 
  // fTopBundles[5] = new CurvedBundle     ( "BundleTopLI_s"        , -124 , 3 , fBundleStartZ , -2  , 25 , -1952 , 19 , 15 , -65 , kFALSE ); 
  // fTopBundles[7] = new CurvedBundle     ( "BundleTopLI_l"        , -124 , 1 , fBundleStartZ , -18 , 25 , -1952 , 19 , 10 , -60 , kFALSE ); 

  // Coordinates of mergin points: 

  // -127 6 // C-Side-Bundle-Top-R-Outer(Back)
  // -124 6 // C-Side-Bundle-Top-R-Inner(Back)
  // -127 2 // C-Side-Bundle-Top-L-Outer(Back)
  // -124 2 // C-Side-Bundle-Top-L-Inner(Back)


}
//_____________________________________________________________________________
void CreateCurvedBundles()
{
  const Int_t fNBundles = 8;
  CurvedBundle * fTopBundles[fNBundles] = {0};
  Double_t fBundleStartZ = -1897;
  // Front and back middle (x>>0) fiber bundles:
  fTopBundles[0] = new CurvedBundle     ( "BundleTop1O"        , -127 , 7 , fBundleStartZ , 18  , 25 , -1956 , 15 , 22 , -75 , kFALSE ); 
  fTopBundles[1] = new CurvedBundle     ( "BundleTop1I"        , -124 , 7 , fBundleStartZ , 18  , 25 , -1952 , 14 , 22 , -75 , kFALSE ); 
  // Front and back middle ( x >0) fiber bundles:
  fTopBundles[2] = new CurvedBundle     ( "BundleTop2O"        , -127 , 5 , fBundleStartZ , 2   , 25 , -1956 , 15 , 18 , -70 , kFALSE ); 
  fTopBundles[3] = new CurvedBundle     ( "BundleTop2I"        , -124 , 5 , fBundleStartZ , 2   , 25 , -1952 , 14 , 18 , -70 , kFALSE ); 
  // Front and back middle ( x <0) fiber bundles:
  fTopBundles[4] = new CurvedBundle     ( "BundleTop3O"        , -127 , 3 , fBundleStartZ , -2  , 25 , -1956 , 20 , 15 , -65 , kFALSE ); 
  fTopBundles[5] = new CurvedBundle     ( "BundleTop3I"        , -124 , 3 , fBundleStartZ , -2  , 25 , -1952 , 19 , 15 , -65 , kFALSE ); 
  // Front and back middle ( x<<0) fiber bundles:
  fTopBundles[6] = new CurvedBundle     ( "BundleTop4O"        , -127 , 1 , fBundleStartZ , -18 , 25 , -1956 , 20 , 10 , -60 , kFALSE ); 
  fTopBundles[7] = new CurvedBundle     ( "BundleTop4I"        , -124 , 1 , fBundleStartZ , -18 , 25 , -1952 , 19 , 10 , -60 , kFALSE ); 
  //

  CurvedBundle * fBottomBundles[fNBundles] = {0};
  // Front and back middle ( x>>0) fiber bundles:
  fBottomBundles[0] = new CurvedBundle  ( "BundleBottom1O"     , -127 , 7 , fBundleStartZ , 18  , 25 , -1956 , 15 , 22 , -75 , kTRUE  ); 
  fBottomBundles[1] = new CurvedBundle  ( "BundleBottom1I"     , -124 , 7 , fBundleStartZ , 18  , 25 , -1952 , 14 , 22 , -75 , kTRUE  ); 
  // Front and back middle ( x >0) fiber bundles:
  fBottomBundles[2] = new CurvedBundle  ( "BundleBottom2O"     , -127 , 5 , fBundleStartZ , 2   , 25 , -1956 , 15 , 18 , -70 , kTRUE  ); 
  fBottomBundles[3] = new CurvedBundle  ( "BundleBottom2I"     , -124 , 5 , fBundleStartZ , 2   , 25 , -1952 , 14 , 18 , -70 , kTRUE  ); 
  // Front and back middle ( x <0) fiber bundles:
  fBottomBundles[4] = new CurvedBundle  ( "BundleBottom3O"     , -127 , 3 , fBundleStartZ , -2  , 25 , -1956 , 20 , 15 , -65 , kTRUE  ); 
  fBottomBundles[5] = new CurvedBundle  ( "BundleBottom3I"     , -124 , 3 , fBundleStartZ , -2  , 25 , -1952 , 19 , 15 , -65 , kTRUE  ); 
  // Front and back middle ( x<<0) fiber bundles:
  fBottomBundles[6] = new CurvedBundle  ( "BundleBottom4O"     , -127 , 1 , fBundleStartZ , -18 , 25 , -1956 , 20 , 10 , -60 , kTRUE  ); 
  fBottomBundles[7] = new CurvedBundle  ( "BundleBottom4I"     , -124 , 1 , fBundleStartZ , -18 , 25 , -1952 , 19 , 10 , -60 , kTRUE  ); 
  //
  // fBottomBundles[0] = new CurvedBundle ("BundleBottom1", -124,5,-1920, 20,25,-1955, 25, 15, -60, kTRUE);
  // fBottomBundles[1] = new CurvedBundle ("BundleBottom2", -124,1,-1920,  2,25,-1955, 25, 15, -60, kTRUE);
  // fBottomBundles[2] = new CurvedBundle ("BundleBottom3", -127,5,-1920, -2,25,-1955, 25, 15, -60, kTRUE);
  // fBottomBundles[3] = new CurvedBundle ("BundleBottom4", -127,1,-1920,-20,25,-1955, 25, 10, -60, kTRUE);

  for (Int_t i=0; i<fNBundles; i++) 
  {
    if (!fTopBundles[i]) continue;
    alice -> AddNode( fTopBundles[i]->GetVolume(), 1);
  }
  for (Int_t i=0; i<fNBundles; i++) 
  {
    if (!fBottomBundles[i]) continue;
    alice -> AddNode( fBottomBundles[i]->GetVolume(), 1);
  }
  // Print Lengths
  for (Int_t i=0; i<fNBundles; i++) 
  {
    if (!fTopBundles[i]) continue;
    Double_t L = fTopBundles[i]->GetTotalLength();
    printf("fTopBundles[%d]: L: %6.3f | 250-L: %6.3f\n", i, L, 250. - L);
  }
  for (Int_t i=0; i<fNBundles; i++) 
  {
    if (!fBottomBundles[i]) continue;
    Double_t L = fBottomBundles[i]->GetTotalLength();
    printf("fBottomBundles[%d]: L: %6.3f | 250-L: %6.3f\n", i, L, 250. - L);
  }
}
//_____________________________________________________________________________

TGeoVolumeAssembly * CreateAD()
{
  printf("AliADv1::CreateAD(): ver=[May 19th, 2017]; contact=[ecalvovi@cern.ch]\n");
  // printf("AD-SHIELDING: %s %d: dy: %f dz: %f off: %f (%p->%s)\n", __FUNCTION__, __LINE__, pmshield_y, pmshield_z, pmshield_offset, fShieldMat, fShieldMat);
  TGeoVolumeAssembly * ad = new TGeoVolumeAssembly("AD");
  Int_t nvertices;
  //
  // Define Rotations used
  //
  Rx90m = new TGeoRotation("Rx90m",   0., -90.,   0.) ;
  Rx90  = new TGeoRotation("Rx90" ,   0.,  90.,   0.) ;
  Rx180 = new TGeoRotation("Rx180",   0., 180.,   0.) ;  //   4    |   1
  Rz180 = new TGeoRotation("Rz180", 180.,   0.,   0.) ;  // --------------->  x  
  Ry180 = new TGeoRotation("Ry180", 180., 180.,   0.) ;  //   3    |   2
  Ry90m = new TGeoRotation("Ry90m",  90., -90., -90.) ;
  Ry90  = new TGeoRotation("Ry90" ,  90.,  90., -90.) ;
  // Get Mediums needed. 
  // TGeoMedium * kMedAlu       = gGeoManager->GetMedium("AD_Alum");   // Aluminium 
  // TGeoMedium * kMedSteelSh   = gGeoManager->GetMedium("AD_ST_C0");  // Stainless Steel 
  // TGeoMedium * kMedVacuum    = gGeoManager->GetMedium("AD_VA_C0");  // Stainless Steel 
  // TGeoMedium * kMedCopper    = gGeoManager->GetMedium("AD_Cu_C0");  // Stainless Steel 
  // TGeoMedium * kMedPVC       = gGeoManager->GetMedium("AD_PVC"  );  // Stainless Steel
  // TGeoMedium * kSensitiveVac = gGeoManager->GetMedium("AD_VA_SENSITIVE");  // Stainless Steel 
  /////////////////////////////////////////////////////////////////////////////
  //
  // Add Cables for A-Side
  // No LHC drawing found!
  //
  new TGeoBBox("shADACablingVBar1", 2.0, 100.0, 2.0);
  new TGeoBBox("shADACablingVBar0", 1.8, 102.0, 1.8);
  TGeoVolume* voADACablingVBar = new TGeoVolume("voADACablingVBar", 
    new TGeoCompositeShape("shADACablingVBar", "shADACablingVBar1-shADACablingVBar0"), kMedSteelSh);
  new TGeoBBox("shADACablingHBar1",  70.0/2.0, 2.0, 2.0);
  new TGeoBBox("shADACablingHBar0",  72.0/2.0, 1.8, 1.8);
  TGeoVolume* voADACablingHBar = new TGeoVolume("voADACablingHBar", 
    new TGeoCompositeShape("shADACablingHBar", "shADACablingHBar1-shADACablingHBar0"), kMedSteelSh);
  voADACablingVBar->SetLineColor(kGray+1);
  voADACablingHBar->SetLineColor(kGray+1);

  /////////////////////////////////////////////////////////////////////////////
  // 
  // Vertical Cables A-Side
  // No LHC drawing found!
  // There are ~ 37 cables in a dx space of 74 cm. That makes cables of 2 cm diameter
  // No info about the cable. Taking as an aproximation a cable from Riyadh catalog:
  // model: 000101xx15 
  // nominal cross section: 1 x 120 
  // number of wires in conductor: 37 
  // diameter of conductor: 14.21 
  // Insulation thickness: 1.6 
  // Overall diameter: 18.8
  //
  // It has 37 17.64
  //
  Double_t rCable     = 2.0/2.;
  Double_t rCableCore = 1.4/2.; 
  Double_t lCable     = 200.;
  // TGeoTube   * shADCableExt  = new TGeoTube   ("shADCableExt" , rCableCore, rCable, lCable/2.);
  // TGeoTube   * shADCableCore = new TGeoTube   ("shADCableCore",         0., rCable, lCable/2.);
  TGeoVolume * voADCableExt  = new TGeoVolume ( "voADCableExt"  , new TGeoTube ( "shADCableExt"  , rCableCore+0.01 , rCable          , lCable/2. )  , kMedPVC) ;
  TGeoVolume * voADCableCore = new TGeoVolume ( "voADCableCore" , new TGeoTube ( "shADCableCore" , 0.              , rCableCore-0.01 , lCable/2. )  , kMedCopper);
  voADCableExt  -> SetLineColor(kGray+3);
  voADCableCore -> SetLineColor(kOrange+2);
  TGeoVolumeAssembly * voADAMagnetCable = new TGeoVolumeAssembly("voADAMagnetCable");
  voADAMagnetCable->AddNode(voADCableExt , 1, Rx90);
  voADAMagnetCable->AddNode(voADCableCore, 1, Rx90);



  TGeoVolumeAssembly * voADAMagnetCableArray= new TGeoVolumeAssembly("voADAMagnetCableArray");
  
  for (Int_t i=0; i<37; i++)
  {
    voADAMagnetCableArray -> AddNode(voADAMagnetCable, i+1, new TGeoTranslation(+rCable + 2.*rCable*i,0,0));
  }
  // Make I-Beam support structure in RB24/3 (A-Side)
  // LHC Drawing: LHCVC2U_0034
  // Dimensions are approximated.
  // Probably a type "HE 100 A" Wide Flange Beam.
  // flange web      :  0.5 cm
  // flange width    : 10.0 cm
  // flange height   :  9.6 cm
  // flange thickness:  0.8 cm
  TGeoXtru * shADsuppIBeam = new TGeoXtru(2);
  shADsuppIBeam->SetNameTitle("shADsuppIBeam","shADsuppIBeam");
  Double_t HIBeam, hIBeam, h2, W, w;
  HIBeam  = 9.6/2.;
  hIBeam  = 0.8;
  h2 = hIBeam;
  w  =  0.5/2.;
  W  = 10.5/2.;
  Double_t xIBeam[] = { -W      , -W            , -w        , -w        , -W            , -W     , W      , W             , w         , w         , W             , W       }; 
  Double_t yIBeam[] = { -HIBeam , hIBeam-HIBeam , h2-HIBeam , HIBeam-h2 , HIBeam-hIBeam , HIBeam , HIBeam , HIBeam-hIBeam , HIBeam-h2 , h2-HIBeam , hIBeam-HIBeam , -HIBeam }; 
  nvertices = sizeof(xIBeam)/sizeof(Double_t);
  shADsuppIBeam->DefinePolygon(nvertices,xIBeam,yIBeam);
  shADsuppIBeam->DefineSection(0, -430., 0., -HIBeam, 1.0); // index, Z position, offset (x,y) and scale for first section
  shADsuppIBeam->DefineSection(1,    0., 0., -HIBeam, 1.0); // idem, second section
  Double_t fOriginBeamCut1[3] = { 3.25, -HIBeam, -4.5};
  Double_t fOriginBeamCut2[3] = {-3.25, -HIBeam, -4.5};
  new TGeoBBox("shADBeamCut1", W/2., HIBeam + hIBeam, 3., fOriginBeamCut1);
  new TGeoBBox("shADBeamCut2", W/2., HIBeam + hIBeam, 3., fOriginBeamCut2);
  TGeoCompositeShape * shADsuppIBeamCutted = new TGeoCompositeShape("shADsuppIBeamCutted", "(shADsuppIBeam-shADBeamCut1)-shADBeamCut2");
  TGeoVolume * voADsuppIBeam = new TGeoVolume("voADsuppIBeam", shADsuppIBeamCutted, kMedSteelSh);

  // Vertical I-Beam:
  TGeoVolume * voADsuppIBeamV = MakeVolIBeam("voADsuppIBeamV"    , kMedSteelSh, 10.5, 9.6, 0.5, 1.6,81.8-9.6);

  /////////////////////////////////////////////////////////////////////////////
  //
  // Make Top Bracket, part of 
  // Sliding support (structure 18 in RB24/3)
  // LHC Drawing: LHCVC2U_0026 (item 5)
  // 
  new TGeoBBox("shADsupp1Box",                  11./2., 15./2., 1. );
  new TGeoBBox("shADsupp1HoleSqr",               8./2.,  4./2., 1.2);  // larger y (4 instead of 3) just in case
  new TGeoTube("shADsupp1HoleCircle",               0.,    0.5, 1.2);
  new TGeoBBox("shADsupp1HoleRectangle",            1.,    0.5, 1.2);
  TGeoArb8 * shADsupp1HoleTrg = new TGeoArb8("shADsupp1HoleTrg", 1.2);
  Float_t trgX[] = { -2.5, 0.0, 2.5,  0.0 };
  Float_t trgY[] = {  0.0, 2.5, 0.0, -2.5 };
  for (Int_t i=0; i<4; i++) 
  {
    shADsupp1HoleTrg -> SetVertex(i,   trgX[i], trgY[i]);
    shADsupp1HoleTrg -> SetVertex(i+4, trgX[i], trgY[i]);
  }
  (new TGeoTranslation("trADsupp1Box",              0., 15./2., 0. ))->RegisterYourself();
  (new TGeoTranslation("trADsupp1HoleSqr",          0.,    14., 0. ))->RegisterYourself();
  (new TGeoTranslation("trADsupp1HoleRndRecL1",     3.,    2.5, 0. ))->RegisterYourself();
  (new TGeoTranslation("trADsupp1HoleRndRecL2",     3.,    5.5, 0. ))->RegisterYourself();
  (new TGeoTranslation("trADsupp1HoleRndRecL3",     3.,    8.5, 0. ))->RegisterYourself();
  (new TGeoTranslation("trADsupp1HoleRndRecR1",    -3.,    2.5, 0. ))->RegisterYourself();
  (new TGeoTranslation("trADsupp1HoleRndRecR2",    -3.,    5.5, 0. ))->RegisterYourself();
  (new TGeoTranslation("trADsupp1HoleRndRecR3",    -3.,    8.5, 0. ))->RegisterYourself();
  (new TGeoTranslation("trADsupp1HoleCircleL",      1.,     0., 0. ))->RegisterYourself();
  (new TGeoTranslation("trADsupp1HoleCircleR",     -1.,     0., 0. ))->RegisterYourself();
  new TGeoCompositeShape("shADsupp1HoleRoundedRectangle", 
      "shADsupp1HoleCircle:trADsupp1HoleCircleL + "
      "shADsupp1HoleCircle:trADsupp1HoleCircleR + "
      "shADsupp1HoleRectangle");
  TGeoCompositeShape * shADsuppTopBracket = new TGeoCompositeShape("shADsuppTopBracket", 
      "shADsupp1Box:trADsupp1Box - "
      "shADsupp1HoleRoundedRectangle:trADsupp1HoleRndRecL1 - "
      "shADsupp1HoleRoundedRectangle:trADsupp1HoleRndRecL2 - "
      "shADsupp1HoleRoundedRectangle:trADsupp1HoleRndRecL3 - "
      "shADsupp1HoleRoundedRectangle:trADsupp1HoleRndRecR1 - "
      "shADsupp1HoleRoundedRectangle:trADsupp1HoleRndRecR2 - "
      "shADsupp1HoleRoundedRectangle:trADsupp1HoleRndRecR3 - "
      "shADsupp1HoleSqr:trADsupp1HoleSqr - "
      "shADsupp1HoleTrg");
  TGeoVolume * voADsuppTopBracket = new TGeoVolume("voADsuppTopBracket", shADsuppTopBracket, kMedAlu);
  

  /////////////////////////////////////////////////////////////////////////////
  //
  // Make support BASE for structure 17 and 18 in RB24/3 (A-Side)
  // LHC Drawing: LHCVC2U_0024
  
  new TGeoBBox("shADsupp17BaseBox",   16./2., 0.65/2., 33./2. );
  new TGeoBBox("shADsupp17VertXBox",  16./2.,   8./2., 0.65/2. );
  new TGeoBBox("shADsupp17VertZBox",  1.3/2.,   8./2.,  25./2. );
  (new TGeoTranslation("trADsupp17BaseBox",  0., 0.65/2., -12.15 ))->RegisterYourself();
  (new TGeoTranslation("trADsupp17VertXBox", 0.,  8.0/2., -12.15 +  8.5 + 0.65/2. ))->RegisterYourself();
  (new TGeoTranslation("trADsupp17VertZBox", 0.,  8.0/2., -12.15 -4. ))->RegisterYourself();
  TGeoCompositeShape * shADsupp17Base = new TGeoCompositeShape("shADsupp17Base", 
      "shADsupp17BaseBox:trADsupp17BaseBox + "
      "shADsupp17VertXBox:trADsupp17VertXBox + "
      "shADsupp17VertZBox:trADsupp17VertZBox");
  TGeoVolume * voADsupp17Base = new TGeoVolume("voADsupp17Base", shADsupp17Base, kMedAlu);
  /////////////////////////////////////////////////////////////////////////////
  //
  // Make support structure 17 in RB24/3 (A-Side)
  // Fixed support
  // LHC Drawing: LHCVC2U_0028
  // LHC Drawing: LHCVC2U_0032
  
  // Main Prop
  TGeoXtru * shADsupp17MainPropLat = new TGeoXtru(2);
  shADsupp17MainPropLat->SetNameTitle("shADsupp17MainPropLat","shADsupp17MainPropLat");
  Double_t x17[] = {  0.0 , 3 , 3.  , 3.   , 3.6  , 5.   , 5.  , 0. };
  Double_t y17[] = {  0.0 , 0 , 0.8 , 7.55 , 7.55 , 8.95 , 28. , 28.};
  nvertices = sizeof(x17)/sizeof(Double_t);
  shADsupp17MainPropLat->DefinePolygon(nvertices,x17,y17);
  shADsupp17MainPropLat->DefineSection(0, 0.00, 0., 0., 1.0); // Z position, offset and scale for first section
  shADsupp17MainPropLat->DefineSection(1, 0.65, 0., 0., 1.0); // idem, second section
  (new TGeoCombiTrans("ctADsupp17BaseLatL"      , 5. - 0.65 , 0.65          , 0. , Ry90 ))->RegisterYourself();
  (new TGeoCombiTrans("ctADsupp17BaseLatR"      , -5.       , 0.65          , 0. , Ry90 ))->RegisterYourself();
  (new TGeoTranslation("trADsupp17MainPropVBox" , 0.        , 28./2. + 0.65 , -0.65/2.  ))->RegisterYourself();
  new TGeoBBox("shADsupp17MainPropVBox", 10./2., 28./2., 0.65/2.);
  TGeoCompositeShape * shADsupp17MainProp = new TGeoCompositeShape("shADsupp17MainProp",
      "shADsupp17MainPropLat:ctADsupp17BaseLatL + "
      "shADsupp17MainPropLat:ctADsupp17BaseLatR + "
      "shADsupp17MainPropVBox:trADsupp17MainPropVBox");
  TGeoVolume * voADsupp17MainProp = new TGeoVolume("voADsupp17MainProp", shADsupp17MainProp, kMedAlu);
  TGeoVolumeAssembly * voADsupp17 = new TGeoVolumeAssembly("voADsupp17");
  voADsupp17 -> AddNode(voADsupp17MainProp, 1);
  voADsupp17 -> AddNode(voADsupp17Base,     1);
  voADsupp17 -> AddNode(voADsuppTopBracket, 1, new TGeoTranslation(0., 15.90, 0.1 + 2./2.));
  


  /////////////////////////////////////////////////////////////////////////////
  //
  // Make support structure 2: part of 
  // Sliding support
  // LHC Drawing: LHCVC2U_0026
  // 
  Float_t htrp   = 4.3;  // real value??
  Float_t trpX[] = {0.  , 5.0 , 4.3,0.  };
  Float_t trpY[] = {htrp, htrp, 0. ,0.  };
  TGeoArb8 * shADsupp2trp = new TGeoArb8("shADsupp2trp", 0.6/2.);
  for (Int_t i=0; i<4; i++) 
  {
    shADsupp2trp -> SetVertex(i,   trpX[i], trpY[i]);
    shADsupp2trp -> SetVertex(i+4, trpX[i], trpY[i]);
  }
  new TGeoBBox("shADsupp2latbox", 5.0/2., (19.-htrp)/2.,  0.6/2.);
  new TGeoBBox("shADsupp2box",    0.6/2.,        19./2., 10.0/2.);
  (new TGeoTranslation("trADsupp2B"   , 0.6/2.,              19./2.,              0.))->RegisterYourself();
  (new TGeoTranslation("trADsupp2L"   , 5.0/2., (19-htrp)/2. + htrp, -10.0/2.+0.6/2.))->RegisterYourself();
  (new TGeoTranslation("trADsupp2R"   , 5.0/2., (19-htrp)/2. + htrp,  10.0/2.-0.6/2.))->RegisterYourself();
  (new TGeoTranslation("trADsupp2trpL",     0.,                  0., -10.0/2.+0.6/2.))->RegisterYourself();
  (new TGeoTranslation("trADsupp2trpR",     0.,                  0.,  10.0/2.-0.6/2.))->RegisterYourself();
  TGeoCompositeShape * shADsupp18MainProp = new TGeoCompositeShape("shADsupp18MainProp", 
      "shADsupp2trp:trADsupp2trpL + "
      "shADsupp2trp:trADsupp2trpR + "
      "shADsupp2latbox:trADsupp2L + "
      "shADsupp2latbox:trADsupp2R + "
      "shADsupp2box:trADsupp2B");
  TGeoVolume * voADsupp18MainProp = new TGeoVolume("voADsupp18MainProp", shADsupp18MainProp, kMedAlu);
  TGeoVolumeAssembly * voADsupp18 = new TGeoVolumeAssembly("voADsupp18");
  voADsupp18 -> AddNode(voADsupp18MainProp, 1, new TGeoCombiTrans (0., 9.65, 0., Ry90m));
  voADsupp18 -> AddNode(voADsupp17Base,     1, Ry180);
  voADsupp18 -> AddNode(voADsuppTopBracket, 1, new TGeoTranslation(0., 16.25, -0.1 - 2./2.));
  /////////////////////////////////////////////////////////////////////////////
  // Define Ion Pump ??
  // Drawing LHCVBU__0052
  // Vacuum - Bellows - U type
  // BODY 1 PORTS
  //
  (new TGeoCombiTrans("ctPumpVB2", 0., -13./2., 6.8-11.5/2., Rx90))->RegisterYourself();
  new TGeoTube("shIonPumpVB1o",  0.0, 10.3 /2., 11.5/2.   );
  new TGeoTube("shIonPumpVB2o",  0.0,  7.0 /2., 13.0/2.   );
  new TGeoTube("shIonPumpVB1i",  0.0, 10.0 /2., 11.5/2.+2.);
  new TGeoTube("shIonPumpVB2i",  0.0,  6.7 /2., 13.0/2.+2.);
  new TGeoCompositeShape("shIonPumpVBo", "shIonPumpVB1o+shIonPumpVB2o:ctPumpVB2");
  //
  // Continue definition of LHCVBU__0052
  //
  new TGeoCompositeShape("shIonPumpVBi", "shIonPumpVB1i+shIonPumpVB2i:ctPumpVB2");
  TGeoShape * sh3 = new TGeoCompositeShape("shIonPumpVB",  "shIonPumpVBo-shIonPumpVBi");
  TGeoVolume * voIonPumpVB = new TGeoVolume("voIonPumpVB", sh3, kMedSteelSh);
  // Variables 
  Double_t alpha, beta, tga2, tga, sa, ca, ctgb, d, Ro, Ri, phi1, dphi, H, L, z;
  // Drawing: LHCVSR__0054
  // Vacuum Screen - RF
  // transition flange
  alpha = 15. * TMath::DegToRad();
  beta  = 15. * TMath::DegToRad();
  tga2  = TMath::Tan(alpha*0.5);
  tga   = TMath::Tan(alpha);
  ctgb  = 1./TMath::Tan(beta );
  ca    = TMath::Cos(alpha    );
  d     = 0.3; // thickness of transition flange
  Double_t h = d/ca; // vertical distance between parallel surfaces tilted alpha degrees
  TGeoPcon * shVSRflange = new TGeoPcon("shVSRflange", 0.0, 360.0, 7);
  Ri = 9.71/2.; Ro = 11.16/2.;
  shVSRflange->DefineSection(0,         0.  ,   Ri, Ro);
  Ri = Ri - 0.25 * tga;
  shVSRflange->DefineSection(1,         0.25,   Ri, Ro);
  Ro = Ri + d*tga;
  shVSRflange->DefineSection(2,         0.25,   Ri, Ro);
  Ri = 6.3/2.; Ro = Ri + h;
  z = (9.71/2. - Ri) / tga;
  shVSRflange->DefineSection(3,            z,   Ri, Ro);
  // 
  Double_t   Dtga = 6.6*tga - 0.5*(9.71-6.3)    ;
  Double_t x = (h - 0.11 - Dtga) / (ctgb - tga) ;
  Double_t y = x * ctgb;
  z  = 6.6 - x;
  Ro = Ri + 0.11 + y;
  shVSRflange->DefineSection(4,            z,   Ri, Ro);
  z  = 6.6; 
  Ro = Ri + 0.11;
  shVSRflange->DefineSection(5,            z,   Ri, Ro);
  z  = 7.1; 
  shVSRflange->DefineSection(6,            z,   Ri, Ro);
  // TGeoVolume * voVSRflange = new TGeoVolume("voVSRflange", shVSRflange, kMedAlu);
  TGeoVolume * voVSRflange = new TGeoVolume("voVSRflange", shVSRflange, kMedCopper);
  //
  // Drawing: LHCVSR__0053
  // Vacuum Screen - RF
  // transition tube
  alpha = 10. * TMath::DegToRad();
  tga2  = TMath::Tan(alpha*0.5);
  tga   = TMath::Tan(alpha);
  ca    = TMath::Cos(alpha    );
  d     = 0.3; // thickness of Vacuum Screen RF
  Ro    = 6.7/2.; //
  Ri    = 0.;
  phi1  = 90. - 15.;
  dphi  = 30.;
  TGeoPcon * shVSR0 = new TGeoPcon("shVSR0", phi1, dphi, 6);
  shVSR0->DefineSection(0,         0.  ,   Ro-0.09, Ro);
  shVSR0->DefineSection(1,         0.45,   Ro-0.09, Ro);
  Ri=Ro-d;
  shVSR0->DefineSection(2,         0.45,      Ri, Ro);
  shVSR0->DefineSection(3, 13.37-d*tga2,      Ri, Ro);
  Ro  = Ri + d/ca;
  shVSR0->DefineSection(4,        13.37,      Ri, Ro);
  Ri += 0.63*tga;
  Ro  = Ri + d/ca;
  shVSR0->DefineSection(5,        14.0 ,      Ri, Ro);
  // printf("  Ro: %8.2f\n", Ro);
  // Make holes 
  new TGeoBBox("shHoleBody"    , 0.15, 0.60, 0.3);
  new TGeoTube("shHoleEnd", 0.  , 0.15, 0.3);
  (new TGeoTranslation("trHoleBody", 0., -0.6, 0.))->RegisterYourself();
  (new TGeoTranslation("trHoleEnd" , 0., -1.2, 0.))->RegisterYourself();
  new TGeoCompositeShape("shHole","shHoleEnd + shHoleEnd:trHoleEnd + shHoleBody:trHoleBody");
  // Single hole made. Now define some combitrans to position holes
  z = 1.3; Ro = (6.7 - d)*0.5;
  (new TGeoCombiTrans("ctHole1", 0., Ro , z, Rx90m))->RegisterYourself();
  z = 1.3 + 1*2.5;
  (new TGeoCombiTrans("ctHole2", 0., Ro , z, Rx90m))->RegisterYourself();
  z = 1.3 + 2*2.5;
  (new TGeoCombiTrans("ctHole3", 0., Ro , z, Rx90m))->RegisterYourself();
  z = 1.3 + 3*2.5;
  (new TGeoCombiTrans("ctHole4", 0., Ro , z, Rx90m))->RegisterYourself();
  // Now make a sector of RF transition tube
  new TGeoCompositeShape("shVSRsec",
   "shVSR0 - (shHole:ctHole1 + shHole:ctHole2 + shHole:ctHole3 + shHole:ctHole4)");
  // Now define rotations for each sector
  TString strSh = "shVSRsec ";
  for (Int_t i=1; i<=11; i++) {
   (new TGeoRotation(Form("rSec%d",i), 30. * i, 0. , 0.))->RegisterYourself();
   strSh+=Form("+ shVSRsec:rSec%d",i);
  }
  // printf("%s\n", strSh.Data());
  TGeoCompositeShape * shVSR = new TGeoCompositeShape("shVSR", strSh.Data());
  // Now assembly the sector to form VSR RF transition tube !
  // TGeoVolume * voVSR = new TGeoVolume("voVSR", shVSR, kMedAlu);
  TGeoVolume * voVSR = new TGeoVolume("voVSR", shVSR, kMedCopper);
  // 
  // Drawing: LHCVSR__0057
  // RF CONTACT
  // 
  Ro = 0.5 * 6.3;
  d  = 0.03;
  Ri = Ro - d;
  // alpha = TMath::ArcSin((7.35-1.75)/Ri);  <-- No!
  H = 0.5 * 6.9 - Ri;
  L = 28. - 7.1 + 0.45 -14. - 1.75;//7.35 - 1.75;
  // Double_t Delta = TMath::Sqrt( L*L + 4.*(H-d)*H );
  Double_t R = TMath::Sqrt((H-d)*(H-d) + L*L);
  alpha = TMath::ASin(d/R) + TMath::ASin((H-d)/R);
  //printf("alpha: %8.2f \n", alpha * TMath::RadToDeg());
  sa = TMath::Sin(alpha);
  ca = TMath::Cos(alpha);
  x = d*sa;
  y = d*ca;
  Double_t R0 =  1.75;
  Double_t R1 = 10.48;
  Double_t R2 =  0.81 + 0.28;
  phi1 = 0.; dphi = 360.; 

  TGeoPcon * shVSRcontact = new TGeoPcon("shVSRcont", phi1, dphi, 6);
  z = 0.;
  shVSRcontact->DefineSection(0, -z, Ri, Ro);
  z = R0;
  shVSRcontact->DefineSection(1, -z, Ri, Ro);
  z  += x;
  Ri += y;
  Ro  = Ri + d;
  shVSRcontact->DefineSection(2, -z, Ri, Ro);
  z  += (R1 - R0) * ca;
  Ri += (R1 - R0) * sa;
  Ro  = Ri + d;
  shVSRcontact->DefineSection(3, -z, Ri, Ro);
  // Last sections (R2)
  Double_t ab = alpha + 21. * TMath::DegToRad();
  Double_t sab = TMath::Sin(ab);
  Double_t cab = TMath::Cos(ab);
  x   = d * sab;
  y   = d/ca - d*cab;
  z  += x;
  Ri += y;
  Ro  = Ri + d/cab;
  shVSRcontact->DefineSection(4, -z, Ri, Ro);
  z  += R2 * cab;
  Ri += R2 * sab;
  Ro  = Ri + d/cab;
  shVSRcontact->DefineSection(5, -z, Ri, Ro);
  TGeoVolume * voVSRcontD = new TGeoVolume("voVSRcontD", shVSRcontact, kMedCopper);
  // TGeoVolume * voVSRcontD = new TGeoVolume("voVSRcontD", shVSRcontact, kMedAlu);
  // Drawing: LHCVSR__0017
  // Vacuum Screen - RF
  // RF Contact flange
  phi1 = 0.;
  dphi = 360.;

  TGeoPcon * shVSRcontFlange = new TGeoPcon("shVSRcontFlange", phi1, dphi, 11);
  Ri = 6.30 /2.; Ro = 11.16 /2.; z = 0;
  shVSRcontFlange->DefineSection( 0, -z, Ri, Ro);
  Ri = 6.30 /2.; Ro = 11.16 /2.; z = 0.1;
  shVSRcontFlange->DefineSection( 1, -z, Ri, Ro);
  Ri = 6.36 /2.; Ro = 11.16 /2.; z = 0.1;
  shVSRcontFlange->DefineSection( 2, -z, Ri, Ro);
  Ri = 6.36 /2.; Ro = 11.16 /2.; z = 0.25;
  shVSRcontFlange->DefineSection( 3, -z, Ri, Ro);
  Ri = 6.36 /2.; Ro =  9.30 /2.; z = 0.25;
  shVSRcontFlange->DefineSection( 4, -z, Ri, Ro);
  Ri = 6.36 /2.; Ro =  9.30 /2.; z = 0.30;
  shVSRcontFlange->DefineSection( 5, -z, Ri, Ro);
  Ri = 6.36 /2.; Ro =  6.90 /2.; z = 0.30;
  shVSRcontFlange->DefineSection( 6, -z, Ri, Ro);
  Ri = 6.36 /2.; Ro =  6.90 /2.; z = 1.10;
  shVSRcontFlange->DefineSection( 7, -z, Ri, Ro);
  Ri = 6.36 /2.; Ro =  6.85 /2.; z = 1.15;
  shVSRcontFlange->DefineSection( 8, -z, Ri, Ro);
  Ri = 6.36 /2.; Ro =  6.50 /2.; z = 1.15;
  shVSRcontFlange->DefineSection( 9, -z, Ri, Ro);
  Ri = 6.36 /2.; Ro =  6.50 /2.; z = 1.30;
  shVSRcontFlange->DefineSection(10, -z, Ri, Ro);
  TGeoVolume * voVSRcontF = new TGeoVolume("voVSRcontF", shVSRcontFlange, kMedCopper);
  // TGeoVolume * voVSRcontF = new TGeoVolume("voVSRcontF", shVSRcontFlange, kMedAlu);
  // Drawing: LHCVBU__0002
  // Bellows + End Parts
  // Vacuum - Bellows - U type

  // First make end part
  phi1 = 0. ; dphi = 360. ;
  TGeoPcon * shVBUend = new TGeoPcon("shVBUend", phi1, dphi, 6);
  Ri = 5.176; Ro = 5.4; z = 0;
  shVBUend->DefineSection( 0, z, Ri, Ro);
  Double_t dz  = 0.03;
  Ri -= dz * TMath::Tan(45. * TMath::DegToRad());
  z  += dz;
  shVBUend->DefineSection( 1, z, Ri, Ro);
  dz  = 0.04;
  Ro -= dz * TMath::Tan(15. * TMath::DegToRad());
  Ri -= dz * TMath::Tan(45. * TMath::DegToRad());
  z  += dz;
  shVBUend->DefineSection( 2, z, Ri, Ro);
  Ro  = 5.250 ;
  Ri  = 5.073;
  z   = 0.103;
  shVBUend->DefineSection( 3, z, Ri, Ro);
  Ro  = 5.25;
  Ri  = 5.02;
  z   = 0.18;
  shVBUend->DefineSection( 4, z, Ri, Ro);
  Ro  = 5.15;
  Ri  = 5.00;
  z   = 0.28;
  shVBUend->DefineSection( 5, z, Ri, Ro);

  TGeoPcon * shVBUtube26mm = new TGeoPcon("shVBUtube26mm", 0., 360., 2);
  shVBUtube26mm->DefineSection( 0, 0.28, 5.0, 5.15);
  shVBUtube26mm->DefineSection( 1, 2.60, 5.0, 5.15);
  TGeoPcon * shVBUtube9mm  = new TGeoPcon("shVBUtube9mm" , 0., 360., 2);
  shVBUtube9mm ->DefineSection( 0, 0.00, 5.0, 5.15);
  shVBUtube9mm ->DefineSection( 1, 0.62, 5.0, 5.15);
  // central part of bellow (TODO: Add plies)
  TGeoPcon * shVBUcent  = new TGeoPcon("shVBUcent" , 0., 360., 6);
  const Int_t nsec = 6;
  Double_t az [nsec] = {0.9, 0.915, 0.915, 11.885, 11.885, 11.9};
  Double_t aRi[nsec] = {5.0, 5.0  , 5.685,  5.685,  5.0  ,  5.0};
  for (Int_t i=0; i<nsec; i++) {
   z=az[i]; Ri = aRi[i]; Ro = 5.7;
   //printf("  i: %2d  z: %8.2f  Ri: %8.2f  Ro: %8.2f\n", i, z, Ri, Ro );
   shVBUcent ->DefineSection( i, z, Ri, Ro);
  }

  ( new TGeoCombiTrans("ctEnd9mm", 0., 0., 0.9, Ry180)) -> RegisterYourself();

  TGeoCompositeShape * shVBU9mm  = new TGeoCompositeShape("shVBU9mm" , "shVBUend:ctEnd9mm + shVBUtube9mm");
  TGeoCompositeShape * shVBU26mm = new TGeoCompositeShape("shVBU26mm", "shVBUend + shVBUtube26mm");

  const Int_t nsec2 = 10;
  Double_t az2 [nsec2] = {0.  , 0.43, 0.43, 0.58, 0.58 , 0.73 , 0.73, 1.05, 1.05, 1.31} ;
  Double_t aRi2[nsec2] = {5.15, 5.15, 5.03, 5.03, 5.455, 5.455, 5.03, 5.03, 5.59, 5.59} ;
  TGeoPcon * shVBUrotFlg  = new TGeoPcon("shVBUrotFlg" , 0., 360., nsec2);
  for (Int_t i=0; i<nsec2; i++) {
   z=az2[i]; Ri = aRi2[i]; Ro = 6.02;
   //printf("  i: %2d  z: %8.2f  Ri: %8.2f  Ro: %8.2f\n", i, z, Ri, Ro );
   shVBUrotFlg ->DefineSection( i, z, Ri, Ro);
  }
  TGeoVolume * voVBUrotFlg = new TGeoVolume("voVBUrotFlg", shVBUrotFlg, kMedSteelSh);
  // TGeoVolume * voVBUrotFlg = new TGeoVolume("voVBUrotFlg", shVBUrotFlg, kMedAlu);

  // Flange 
  TGeoPcon * shVBUflg  = new TGeoPcon("shVBUflg" , 0., 360., 4);
  z  = 0; 
  Ri = 6.02; Ro = 7.6; 
  shVBUflg->DefineSection(0, z, Ri, Ro);
  z  = 1.31; 
  shVBUflg->DefineSection(1, z, Ri, Ro);
  z  = 1.31; 
  Ri = 5.15; 
  shVBUflg->DefineSection(2, z, Ri, Ro);
  z  = 2.0;
  shVBUflg->DefineSection(3, z, Ri, Ro);
  TGeoVolume * voVBUflg = new TGeoVolume("voVBUflg", shVBUflg, kMedSteelSh);
  // TGeoVolume * voVBUflg = new TGeoVolume("voVBUflg", shVBUflg, kMedAlu);

   
  // Add the metal plate at the end of Absorber
  // Plate:  80 cm x 80 cm x 1.95 cm (Thickness is aproximatted) (ernesto.calvo@pucp.edu.pe)
  // The End of the concrete abosorber is at kZendAbs = 1880.75

  new TGeoBBox("shBasePlate", 80./2., 80./2.,  1.95/2.);
  new TGeoTube("shHolePlate",  0.   , 12.3  ,  1.95   );
  TGeoVolume* voSaa3EndPlate  =  new TGeoVolume("voYSaa3EndPlate",
      new TGeoCompositeShape("shYSaa3EndPlate","shBasePlate-shHolePlate"),
      kMedSteelSh);
  //
  // Add Rods
  //
  // dimensions of rods
  const Double_t dzRodL = 27.0;
  const Double_t dzRodA =  4.3;
  const Double_t dzRodB =  1.3;
  new TGeoTube("shLargeRod", 0.,   1.6/2.,  dzRodL/2.);
  new TGeoTube("shRodA",     0.,   3.0/2.,  dzRodA/2.);
  new TGeoTube("shRodB",     0.,   2.3/2.,  dzRodB/2.);
  //
  ( new TGeoTranslation("trRod1", 0., 0., -dzRodL/2. + dzRodA/2.)          )->RegisterYourself();
  ( new TGeoTranslation("trRod2", 0., 0., -dzRodL/2. + dzRodA + dzRodB/2.) )->RegisterYourself();
  ( new TGeoTranslation("trRod3", 0., 0.,  dzRodL/2. - dzRodB/2.)          )->RegisterYourself();
  TGeoVolume * voSaa3Rod = new TGeoVolume("YSAA3Rod",
      new TGeoCompositeShape("shLargeRod+shRodA:trRod1 + shRodB:trRod2 + shRodB:trRod3"),
      kMedSteelSh);
  //
  // Define Valve support (VS)  (ernesto.calvo@pucp.edu.pe) 
  //
  const Double_t dyVS =  5.5; 
  const Double_t dxVS = 30.0; 
  const Double_t dzVS =  1.0;
  TGeoVolume * voVS = new TGeoVolume("voVS",
      new TGeoBBox("shVS", dxVS/2., dyVS/2., dzVS/2.),
      kMedSteelSh);
  voVS->SetLineColor(kGray);

  // Add Valve (Valve is divided in parts VA,VB,VC and VD)  (ernesto.calvo@pucp.edu.pe) 
  TGeoVolumeAssembly * voValve = new TGeoVolumeAssembly("voValve");
  //
  // Define volume VA  (ernesto.calvo@pucp.edu.pe) 
  const Double_t dxVA  = 20.3; 
  const Double_t dyVA  = 48.0; 
  const Double_t dzVA  =  6.0; // Width 
  const Double_t dz2VA =  8.5; // Full width including protuding boxes
  // Valve position  (ernesto.calvo@pucp.edu.pe) 
  const Double_t zPosValve =  kZbegValve + dz2VA/2.;
  //
  new TGeoBBox("shVAbox",       dxVA/2., dyVA/2.,     dzVA/2.);
  new TGeoBBox("shVAHbox",  -1.+dxVA/2.,   3./2.,    dz2VA/2.);
  new TGeoTube("shVAC",              0.,     7.9,    dz2VA/2.);
  new TGeoTube("shVACh",             0.,     5.3,    dz2VA   );
  // translation for shVAHbox (ernesto.calvo@pucp.edu.pe) 
  ( new TGeoTranslation("trVAH1", 0.,  12.75, 0.) )->RegisterYourself();
  ( new TGeoTranslation("trVAH2", 0., -12.75, 0.) )->RegisterYourself();

  TGeoVolume * voValveVA = new TGeoVolume("voValveVA",
      new TGeoCompositeShape("(shVAbox + shVAHbox:trVAH1 + shVAHbox:trVAH2 + shVAC)-shVACh"),
      kMedSteelSh);
  voValve->AddNode(voValveVA, 1, 0);
  // Define Vacuum Hole of Valve
  TGeoTube   * shVACvach   = new TGeoTube("shVACvach", 0., 5.3, dz2VA/2.);
  TGeoVolume * voValveVAvh = new TGeoVolume("voValveVAvacuum", shVACvach, kMedVacuum);
  voValve->AddNode(voValveVAvh,1,0);
  // Also add valve Support (ernesto.calvo@pucp.edu.pe) 
  voValve->AddNode(voVS, 1, new TGeoTranslation(0.,  12.75, -dz2VA/2.-dzVS/2.));
  voValve->AddNode(voVS, 2, new TGeoTranslation(0., -12.75, -dz2VA/2.-dzVS/2.));

  // Define Volume VB (ernesto.calvo@pucp.edu.pe) 
  const Double_t dxVB = 23.5; 
  const Double_t dyVB =  5.0; 
  const Double_t dzVB =  9.4; 
  TGeoVolume * voValveVB = new TGeoVolume("voValveVB", 
      new TGeoBBox("shVBbox", dxVB/2., dyVB/2., dzVB/2.),
      kMedSteelSh);
  voValve->AddNode(voValveVB, 1, new TGeoTranslation(  0., dyVA/2. +dyVB/2. , 0));
  // Define Volume VC (ernesto.calvo@pucp.edu.pe) 
  const Double_t R1VC  =  4.5 /2.;
  const Double_t R2VC  =  8.1 /2.;
  const Double_t dy1VC  = 10.0;
  const Double_t dy2VC =  0.75;
  new TGeoTube("shVC1",      0.,   R1VC, dy1VC/2.);
  new TGeoTube("shVC2",      0.,   R2VC, dy2VC/2.);
  ( new TGeoTranslation("trVC21", 0., 0.,  dy1VC/2. - dy2VC/2.) )->RegisterYourself();
  ( new TGeoTranslation("trVC22", 0., 0., -dy1VC/2. + dy2VC/2.) )->RegisterYourself();
  TGeoVolume * voValveVC = new TGeoVolume("voValveVC",
      new TGeoCompositeShape("shVC1  + shVC2:trVC21 + shVC2:trVC22"),
      kMedSteelSh);
  voValve->AddNode(voValveVC, 1, new TGeoCombiTrans( 
        0., dyVA/2. + dyVB + dy1VC/2. , 0, Rx90) );
  // Define volume VD (ernesto.calvo@pucp.edu.pe) 
  const Double_t dxVD = 15.9;
  const Double_t dyVD = 23.0;
  const Double_t dzVD = 14.0;
  TGeoVolume * voValveVD = new TGeoVolume("voValveVD",
      new TGeoBBox("shVD", dxVD/2., dyVD/2., dzVD/2.),
      kMedSteelSh);
  voValve->AddNode(voValveVD, 1, 
      new TGeoTranslation( 1.25, dyVA/2. + dyVB + dy1VC + dyVD/2. , 0) );

  //
  // Define volume Front Bar (ernesto.calvo@pucp.edu.pe) 
  //
  const Double_t dxF  = 67.4; 
  const Double_t dyF  =  4.0; 
  const Double_t dzF  =  2.0; 
  const Double_t R1FC =  8.1;
  const Double_t R2FC = 11.5;
  const Double_t dxFA  = (dxF-R1FC)/2.; 
 
  new TGeoBBox("shFA", dxFA/2., dyF/2., dzF/2.);
  new TGeoTube("shFC",      0.,   R2FC, dzF/2.);
  new TGeoTube("shFCH",     0.,   R1FC, 2.*dzF);
  ( new TGeoTranslation("trFA1",  R1FC +dxFA/2., 0., 0.) )->RegisterYourself();
  ( new TGeoTranslation("trFA2", -R1FC -dxFA/2., 0., 0.) )->RegisterYourself();
  TGeoVolume * voFrontBar = new TGeoVolume("voFrontBar",
      new TGeoCompositeShape("shFA:trFA1 + shFA:trFA2 + (shFC - shFCH)"),
      kMedSteelSh);
  // Make Lateral Bars
  const Double_t kdzLatBar = 22.9;
  TGeoVolume * voLatBar = new TGeoVolume("voLatBar", 
    new TGeoTube("shLatBar",0., dyF/2., kdzLatBar/2.), kMedSteelSh);
  //
  // Define Compensator Magnet coils
  //
  dz = 12.5;
  Ro = 15.0;
  Ri = 1.9;
  new TGeoTubeSeg("shCoilRo", 0., Ro, dz/2., 90., 180.);
  new TGeoTubeSeg("shCoilRi", 0., Ri, dz   , 90., 185.);
  (new TGeoTranslation("trCoilRo", Ro, 10.4, 0.)) -> RegisterYourself();
  (new TGeoTranslation("trCoilRi", Ro,  9.6, 0.)) -> RegisterYourself();
  new TGeoBBox("shBoxCoilRo", 15.0/2., 10.4/2., dz/2.);
  new TGeoBBox("shBoxCoilRi",  1.9/2.,  9.6   , dz   );
  (new TGeoTranslation("trBoxCoilRo",       15.0/2., 10.4/2., 0.)) -> RegisterYourself();
  (new TGeoTranslation("trBoxCoilRi", 13.1 + 1.9/2.,  0.    , 0.)) -> RegisterYourself();
  new TGeoBBox("shBoxCoil0", 10./2., 30., dz);
  (new TGeoTranslation("trBoxCoil0",  14.6 + 10./2.,  0.    , 0.)) -> RegisterYourself();
  strSh  = "";
  strSh += "(shCoilRo:trCoilRo + shBoxCoilRo:trBoxCoilRo) - ";
  strSh += "(shCoilRi:trCoilRi + shBoxCoilRi:trBoxCoilRi +"  ;
  strSh += " shBoxCoil0:trBoxCoil0 )"  ;
  TGeoCompositeShape * shCoil = new TGeoCompositeShape("shCoil0", strSh);
  TGeoVolume * voCoil = new TGeoVolume("voCoil", shCoil, kMedSteelSh);

  // 
  // ALUMINIUM PLATES 
  //

  // Shape for aluminium Plate separating cavern and LHC tunnel
  const Double_t dAlWallThick = 0.5; // thickness of aluminium plates (cm)
  //
  // RB24/26 Tunnel Floor 
  R   = 220.;
  // h   = 140.;
  // phi = TMath::ACos(h / r);
  // xl  = r * TMath::Sin(phi);
  // dr  = 1600.;
  // dh  = dr * TMath::Cos(phi);
  // dl  = dr * TMath::Sin(phi);

  new TGeoTube("shWallBase",    0.,    R, dAlWallThick*0.5);  // base shape for shWallBigPlate
  new TGeoBBox("shWallCutBot",270., 110., dAlWallThick    );  // to be substracted from base
  Double_t origin[3] = {+120.+1.5+3.25, 0., 0.};
  new TGeoBBox("shWallBundleHole", 3.25,  7.5, dAlWallThick, origin);  // to be substracted from base
  // Translation for cutting circular and square hole in the plates
  (new TGeoTranslation("trAntiBeamAxis",   -70.,               40., 0.)) -> RegisterYourself();
  (new TGeoTranslation("trHUWAT3",         -70., -110. - 140. +40., 0.)) -> RegisterYourself();

  //
  //  Wall Big Aluminium Plate with Squared Hole 
  //
  const Double_t dSqrHoleSide = 33.0; // Side
  new TGeoBBox("shWallSqrHole", dSqrHoleSide*0.5, dSqrHoleSide*0.5, dAlWallThick);
  strSh  = ""; 
  strSh += "shWallBase:trAntiBeamAxis - ";
  strSh += " ( shWallCutBot:trHUWAT3" ;
  strSh += " + shWallSqrHole";
  strSh += " + shWallBundleHole)";
  TGeoVolume* voWallBigPlate = new TGeoVolume("voWallBigPlate", 
    new TGeoCompositeShape("shWallBigPlate", strSh), kMedAlu ); // commented on: 2016-Dec-05
    // new TGeoCompositeShape("shWallBigPlate", strSh), kMedSteelSh); // Added     on: 2016-Dec-05
  //
  // Wall Squared Aluminium Plate 
  //
  const Double_t dCircHoleRad = 9.5; // Radius
  new TGeoTube("shCircHole", 0., dCircHoleRad, dAlWallThick);
  // Make holes for bars
  new TGeoTube("shRodHole", 0.,   1.7/2.,  2*dAlWallThick);
  ( new TGeoTranslation("trWallRod1",  12.5, -12.75, 0.)) -> RegisterYourself();
  ( new TGeoTranslation("trWallRod2",  12.5,  12.75, 0.)) -> RegisterYourself();
  ( new TGeoTranslation("trWallRod3", -12.5, -12.75, 0.)) -> RegisterYourself();
  ( new TGeoTranslation("trWallRod4", -12.5,  12.75, 0.)) -> RegisterYourself();
  new TGeoBBox("shWallSqrPlateBase", dSqrHoleSide*0.5 + 5.0, dSqrHoleSide*0.5 + 5.0, dAlWallThick/2.);
  strSh  = " ((((";
  strSh += " ( shWallSqrPlateBase - shCircHole )"; 
  strSh += " - shRodHole:trWallRod1)"  ;
  strSh += " - shRodHole:trWallRod2)"  ;
  strSh += " - shRodHole:trWallRod3)"  ;
  strSh += " - shRodHole:trWallRod4)"  ;
  TGeoVolume* voWallSqrPlate = new TGeoVolume("HUWAT_AlWall02", 
    new TGeoCompositeShape("shWallSqrPlate", strSh ), kMedAlu);
  // ==========================================================================
  //
  // Define Mother Vacuum volume of VMAOI  (need shIonPumpVBo)
  //
  // ==========================================================================
  const Double_t kdzMoFlange   =  2.0;
  const Double_t kdzMoBellow   = 15.6;
  const Double_t kdzTTube = 11.5; // Bellow starts here
  const Double_t kziTTube  =  1.0; // Ion Pum Tube starts here
  new TGeoTube( "shMoFlange", 0., 15.2/2., kdzMoFlange/2.0);
  new TGeoTube( "shMoBellow", 0., 11.4/2., kdzMoBellow/2.0);
  (new TGeoTranslation("trMoFlange1", 0., 0.,        0.5*kdzMoFlange     )) -> RegisterYourself();
  (new TGeoTranslation("trMoFlange2", 0., 0., 28.0 - 0.5*kdzMoFlange     )) -> RegisterYourself();
  (new TGeoTranslation("trMoBellow" , 0., 0., 28.0 - 0.5*kdzMoBellow     )) -> RegisterYourself();
  (new TGeoTranslation("trMoTTube"  , 0., 0., kziTTube +  0.5*kdzTTube   )) -> RegisterYourself();
  
  TGeoCompositeShape * shMoVMAOI =  
      new TGeoCompositeShape("shMoVMAOI", 
      "shMoFlange:trMoFlange1 + shMoBellow:trMoBellow + shMoFlange:trMoFlange2 + shIonPumpVBo:trMoTTube");
  TGeoVolume * voMoVMAOI = new TGeoVolume("voMoVMAOI", shMoVMAOI, kMedVacuum);
  //
  // TGeoVolume * voVBU9mm  = new TGeoVolume("voVBU9mm" , shVBU9mm , kMedAlu);
  // TGeoVolume * voVBU26mm = new TGeoVolume("voVBU26mm", shVBU26mm, kMedAlu);
  // TGeoVolume * voVBUcent = new TGeoVolume("voVBUcent", shVBUcent, kMedAlu);
  //
  TGeoVolume * voVBU9mm  = new TGeoVolume("voVBU9mm" , shVBU9mm , kMedSteelSh);
  TGeoVolume * voVBU26mm = new TGeoVolume("voVBU26mm", shVBU26mm, kMedSteelSh);
  TGeoVolume * voVBUcent = new TGeoVolume("voVBUcent", shVBUcent, kMedSteelSh);
  voVSR      ->SetLineColor(kViolet+6);
  voIonPumpVB->SetLineColor(kViolet+6);
  voVSRflange->SetLineColor(kViolet+6);
  voVSRcontD ->SetLineColor(kViolet+6);
  voVSRcontF ->SetLineColor(kViolet+6);
  voVBUrotFlg->SetLineColor(kViolet+6); 
  voVBUrotFlg->SetLineColor(kViolet+6); 
  voVBUflg   ->SetLineColor(kViolet+6); 
  voVBU9mm   ->SetLineColor(kViolet+6);
  voVBU26mm  ->SetLineColor(kViolet+6);
  voVBUcent  ->SetLineColor(kViolet+6);
  //         
  voMoVMAOI->AddNode(voVSR      ,1, new TGeoTranslation(0., 0., 7.1 - 0.45));
  voMoVMAOI->AddNode(voIonPumpVB,1, new TGeoTranslation(0., 0., 1 + 11.5/2.));
  voMoVMAOI->AddNode(voVSRflange,1, new TGeoTranslation(0.,0.,0.));
  voMoVMAOI->AddNode(voVSRcontD ,1, new TGeoTranslation(0.,0.,28.));
  voMoVMAOI->AddNode(voVSRcontF ,1, new TGeoTranslation(0.,0.,28.));
  z = 1.0 + 11.5;
  voMoVMAOI->AddNode( voVBU9mm   , 1, new TGeoTranslation(0.,0., z) );
  voMoVMAOI->AddNode( voVBU26mm  , 1, new TGeoTranslation(0.,0., z + 11.9) );
  voMoVMAOI->AddNode( voVBUcent  , 1, new TGeoTranslation(0.,0., z) );
  voMoVMAOI->AddNode( voVBUrotFlg, 1, new TGeoCombiTrans (0.,0.,1.31, Ry180) );
  voMoVMAOI->AddNode( voVBUrotFlg, 2, new TGeoTranslation(0.,0.,28. - 1.31) );
  voMoVMAOI->AddNode( voVBUflg   , 1, new TGeoTranslation(0.,0.,0.) );
  voMoVMAOI->AddNode( voVBUflg   , 2, new TGeoCombiTrans (0.,0.,28., Ry180) );
  // ==========================================================================
  //
  // AD Support structure by Pieter Ijzerman
  // ecalvovi@cern.ch
  // ==========================================================================
  nvertices=0;
  // Cover plate_______________________________________________________________
  TGeoXtru * shADcoverplate = new TGeoXtru(2);
  shADcoverplate->SetNameTitle("shADcoverplate","shADcoverplate");
  Double_t y1[] = {  0.0, 18.50, 18.50, 22.50, 22.50, 18.50, 18.50, 22.50, 22.50, 18.50, 18.50,   .00 ,  .00, 15.25, 15.25,  .00 }; 
  Double_t x1[] = {  0.0,   .00,  5.15,  5.15, 17.15, 17.15, 24.25, 24.25, 36.25, 36.25, 41.40, 41.40 ,35.70, 35.70,  5.70, 5.70 }; 
  nvertices = sizeof(x1)/sizeof(Double_t);
  shADcoverplate->DefinePolygon(nvertices,x1,y1);
  shADcoverplate->DefineSection(0, -0.1, -20.7, 0.0, 1.0); // Z position, offset and scale for first section
  shADcoverplate->DefineSection(1,  0.1, -20.7, 0.0, 1.0); // -''- secons section

  // Horizontal side___________________________________________________________
  TGeoXtru * shADhorizontalside = new TGeoXtru(2);
  shADhorizontalside->SetNameTitle("shADhorizontalside","shADhorizontalside");
  Double_t x2[] = {  0.0,  .00, 4.80, 4.80, 7.20, 7.20, 12.00, 12.00 };
  Double_t y2[] = {  0.0, 5.66, 5.66, 1.16, 1.16, 5.66,  5.66,   .00 };
  nvertices = sizeof(x2)/sizeof(Double_t);
  shADhorizontalside->DefinePolygon(nvertices,x2,y2);
  shADhorizontalside->DefineSection(0, -0.4, -6.0, 0.0, 1.0); // Z position, offset and scale for first section
  shADhorizontalside->DefineSection(1, +0.4, -6.0, 0.0, 1.0); // -''- secons section

  TGeoBBox * shADsidebox = new TGeoBBox("shADsidebox", 0.4, 18.55/2., 5.66/2.);
  TGeoVolume * voADsidebox = new TGeoVolume("voADsidebox", shADsidebox, kMedAlu);


  TGeoVolume * voADhorizontalside = new TGeoVolume("voADhorizontalside", shADhorizontalside, kMedAlu);
  TGeoVolume * voADcoverplate = new TGeoVolume("voADcoverplate", shADcoverplate, kMedAlu);
  //
  TGeoVolume *voADsupport = new TGeoVolumeAssembly("voADsupport"); 
  voADsupport->AddNode(voADcoverplate,  1, new TGeoTranslation( 0., 0., -5.66/2.-0.23));
  voADsupport->AddNode(voADcoverplate,  2, new TGeoTranslation( 0., 0., +5.66/2.+0.23));
  voADsupport->AddNode(voADhorizontalside,  1, new TGeoCombiTrans( -6.0 - 7.1/2., 22.5-0.4, -5.66/2., Rx90));
  voADsupport->AddNode(voADhorizontalside,  2, new TGeoCombiTrans( +6.0 + 7.1/2., 22.5-0.4, -5.66/2., Rx90));
  voADsupport->AddNode(voADsidebox,  1, new TGeoTranslation( -20.7 +0.4, 18.55/2., 0.));
  voADsupport->AddNode(voADsidebox,  2, new TGeoTranslation( +20.7 -0.4, 18.55/2., 0.));
  // Add Color
  voADcoverplate     ->SetLineColor(kGray + 1);
  voADhorizontalside ->SetLineColor(kGray + 1);
  voADsidebox        ->SetLineColor(kGray + 1);
  voADsupport->SetLineColor(kGray+1);

  // ==========================================================================
  //
  // Define ADA
  //
  // ==========================================================================

  
  // Get medium for ADA
  TGeoMedium * medADASci        = gGeoManager->GetMedium("AD_BC404"); // AD Scin.
  // TGeoMedium * medADALG      = gGeoManager->GetMedium("AD_PMMA");  // lightGuide
  // TGeoMedium * medADAPMGlass = gGeoManager->GetMedium("AD_Glass"); // Glass for Aluminium simulation
  // TGeoMedium * medADAPMAlum  = gGeoManager->GetMedium("AD_Alum");
  
  // Get Medium for ADC 
  TGeoMedium * medADCSci     = gGeoManager->GetMedium("AD_BC404");
  // TGeoMedium * medADCLG      = gGeoManager->GetMedium("AD_PMMA");
  // TGeoMedium * medADCPMGlass = gGeoManager->GetMedium("AD_Glass");
  // TGeoMedium * medADCPMAlum  = gGeoManager->GetMedium("AD_Alum");
  
  // ADA Scintillator Pad 
  const Double_t kADACellSideY = 21.6;
  const Double_t kADACellSideX = 18.1;
  // ADC Scintillator Pad 
  const Double_t kADCCellSideY = 21.6;
  const Double_t kADCCellSideX = 18.1;
  // WLS bar          :  0.40 cm ( 4.0 mm )
  // Wrapping         :  0.20 cm ( 2.0 mm )
  // Aluminnized Mylar:  0.01 cm ( 0.1 mm )
  // Fishing line     :  0.04 cm ( 0.4 mm )
  // total shift on X :  0.65 cm
  // total shift on Y :  0.21 cm
  const Double_t kShiftX       =  0.54;
  const Double_t kShiftY       =  0.10;
  const Double_t kADACelldz    =  2.54;
  const Double_t kADCCelldz    =  2.54;
  const Double_t kADABeamPipeR =  6.20; // Radius of beam pipe hole for ADA (Diameter  12.4 cm)
  const Double_t kADCBeamPipeR =  3.70; // Radius of beam pipe hole for ADC (Diameter   7.4 cm)
  const Int_t    kColorADA     = kGreen;
  const Int_t    kColorADC     = kGreen;
  Double_t X = kShiftX + kADACellSideX * 0.5;
  Double_t Y = kShiftY + kADACellSideY * 0.5;
  Double_t WLS_dx =  0.4;
  Double_t WLS_dz =  2.5;
  Double_t WLS_SideA_Long_dy  = 24.20; // 24.2;
  Double_t WLS_SideC_Long_dy  = 24.20; // 24.2;
  Double_t WLS_SideA_Short_dy = 18.20; // 18.41; 
  Double_t WLS_SideC_Short_dy = 20.70; // 20.91; 
  // Creating ADA WLS bars_____________________________________________________
  TGeoVolume * vADA_WLS_s = new TGeoVolume( "ADAWLSshort", 
      new TGeoBBox( "shADAWLSbarShort" , WLS_dx/2.0, WLS_SideA_Short_dy/2.0, WLS_dz/2.0),
      medADASci);      
  TGeoVolume * vADA_WLS_l = new TGeoVolume( "ADAWLSlong" , 
      new TGeoBBox( "shADAWLSbarLong"  , WLS_dx/2.0, WLS_SideA_Long_dy /2.0, WLS_dz/2.0),
      medADASci);      
  vADA_WLS_l->SetLineColor( kRed );
  vADA_WLS_s->SetLineColor( kRed );
  // Creating ADC WLS bars_____________________________________________________
  TGeoVolume * vADC_WLS_s = new TGeoVolume( "ADCWLSshort", 
      new TGeoBBox( "shADCWLSbarShort" , WLS_dx/2.0, WLS_SideC_Short_dy/2.0, WLS_dz/2.0),
      medADCSci);      
  TGeoVolume * vADC_WLS_l = new TGeoVolume( "ADCWLSlong" , 
      new TGeoBBox( "shADCWLSbarLong"  , WLS_dx/2.0, WLS_SideC_Long_dy /2.0, WLS_dz/2.0),
      medADCSci);      
  vADC_WLS_l->SetLineColor(kRed);
  vADC_WLS_s->SetLineColor(kRed);
  // Make ADA scintillator pad_________________________________________________
  new TGeoBBox( "shADAbox" , kADACellSideX/2.0, kADACellSideY/2.0, kADACelldz/2.0 );
  new TGeoTube( "shADAHole",               0. , kADABeamPipeR    , kADACelldz     );
  ( new TGeoTranslation("trADAbox", X, Y, 0.)) -> RegisterYourself();
  // 
  TGeoVolume * vADA1 = new TGeoVolume( "ADApad", 
    new TGeoCompositeShape("shADApad", "shADAbox:trADAbox-shADAHole"), medADASci );      
  vADA1->SetLineColor( kColorADA );
  
  TGeoVolume *secADA  = new TGeoVolumeAssembly( "ADAsec" ); 
  // Add PAD
  Double_t fX_ADA_WLS_s = 0.1 + WLS_dx/2.0;
  Double_t fX_ADA_WLS_l = kShiftX + WLS_dx/2.0 + kADACellSideX + 0.04;
  secADA->AddNode( vADA1, 1, 0); 
  secADA->AddNode( vADA_WLS_s, 1, new TGeoTranslation(fX_ADA_WLS_s, kADABeamPipeR + WLS_SideA_Short_dy/2.0, 0.0) ); 
  secADA->AddNode( vADA_WLS_l, 1, new TGeoTranslation(fX_ADA_WLS_l,        kShiftY + WLS_SideA_Long_dy/2.0, 0.0) ); 

  /// Assembling ADA adding 4 sectors                                       //  Sectors
  TGeoVolume *vADAarray = new TGeoVolumeAssembly( "ADA" );                  //        ^ y
  vADAarray->AddNode( secADA, 1 );                                          //        |   
  vADAarray->AddNode( secADA, 2, Ry180 );                                   //   2    |   1
  vADAarray->AddNode( secADA, 3, Rz180 );                                   // --------------->  x     
  vADAarray->AddNode( secADA, 4, Rx180 );                                   //   3    |   4
  // 
  // PMT-BOX A-Side
  //
  const Float_t  kPosADA = 1696.67;  // z-center of assembly (cm) New, according to Survey by F. Klumb and E.Calvo 2015 Sept 4th.
  const Double_t kzEndPMboxA = 1701.65;
  const Double_t kdzPMboxA   = 15.6;
  Float_t thickbox = 0.3; // in cm
  Float_t pmbox_x = 60.;
  Float_t pmbox_y = 45.;
  Float_t pmbox_z = kdzPMboxA;
  // Float_t pmbox_z = 10.;
  Double_t obox3[3] = {0., -2., 0.};
  new TGeoBBox("shPMTbox1" , pmbox_x/2.             , pmbox_y/2.           , pmbox_z/2.           );
  new TGeoBBox("shPMTbox2" , pmbox_x/2. -  thickbox , pmbox_y/2. -thickbox , pmbox_z/2. -thickbox );
  new TGeoBBox("shPMTbox3" , pmbox_x/2. -2*thickbox , pmbox_y/2. -thickbox , (pmbox_z-4.)/2.         , obox3);
  TGeoVolume * voPMTbox = new TGeoVolume("voPMTbox", new TGeoCompositeShape("shPMTbox","(shPMTbox1-shPMTbox2)-shPMTbox3"), kMedAlu);
  voPMTbox -> SetLineColor(kGray + 1);
  //
  // PMT's  A-Side
  //
  Double_t fPMTdz = 2.0; // Length
  Double_t fPMTDi = 2.7; // Diameter
  TGeoVolume * voADApmt = new TGeoVolume("voADApmt", new TGeoTube("shADApmt", 0., fPMTDi/2., fPMTdz/2.), medADASci);
  voADApmt -> SetLineColor(kGray);
  // PMT boxes 
  // Double_t aX_PMT[8] = {6, -12, -12, 6, 12, -6, -6, 12};
  // Double_t aY_PMT[8] = {
  //    82.0+ fPMTdz/2.,
  //    82.0+ fPMTdz/2.,
  //  -124.6- fPMTdz/2.,
  //  -124.6- fPMTdz/2.,
  //    82.0+ fPMTdz/2.,
  //    82.0+ fPMTdz/2.,
  //  -124.6- fPMTdz/2.,
  //  -124.6- fPMTdz/2.
  // };
  Double_t aX_PMT[8] = { 6.5, -11.5, -12.0, 6.0, 12.5, -6.0, -6.0, 12.0 };
  Double_t aY_PMT[8] = {
      +77.2 + fPMTdz/2. , 
      +79.0 + fPMTdz/2. , 
     -120.5 - fPMTdz/2. , 
     -121.9 - fPMTdz/2. , 
      +78.1 + fPMTdz/2. , 
      +77.4 + fPMTdz/2. , 
     -123.0 - fPMTdz/2. , 
     -121.0 - fPMTdz/2.
  };
 

  // Layer closer to IP is shifted to the right in picture (i.e. towards negative x axis) 
  // sector 1: Inner PMT @ x=  6cm 
  // sector 2: Inner PMT @ x=-12cm 
  // sector 3: Inner PMT @ x=-12cm
  // sector 4: Inner PMT @ x=  6cm 
  //
  // sector 1: Outer PMT @ x= 12cm 
  // sector 2: Outer PMT @ x= -6cm 
  // sector 3: Outer PMT @ x= -6cm
  // sector 4: Outer PMT @ x= 12cm 
  //
  TGeoVolumeAssembly * voADAPMTarray = new TGeoVolumeAssembly("voADAPMTarray");
  // Add PMT Boxes

  Double_t zCen_PMbox_wrt_AD = (kzEndPMboxA - kPosADA) - 0.5 * kdzPMboxA;
  Double_t zBeg_PMbox_wrt_AD = (kzEndPMboxA - kPosADA) - 1.0 * kdzPMboxA;
  Double_t zEnd_PMbox_wrt_AD = (kzEndPMboxA - kPosADA);
  voADAPMTarray->AddNode( voPMTbox,    1, new TGeoTranslation( 4.5,   51.+ pmbox_y/2., zCen_PMbox_wrt_AD)); 
  voADAPMTarray->AddNode( voPMTbox,    2, new TGeoCombiTrans ( 3.5, -145.+ pmbox_y/2., zCen_PMbox_wrt_AD, Rz180)); 


  // PMT's
  for (Int_t i=0; i<8; i++) {
    Double_t Z1 = 0;
    if (i<4) Z1 = +kADACelldz/2. +0.23;
    else     Z1 = +kADACelldz/2. +0.23; 
    voADAPMTarray->AddNode( voADApmt, i+1, new TGeoCombiTrans(aX_PMT[i], aY_PMT[i], Z1, Rx90)); 
  }
  //
  // Add the Fiber-Bundles (A-Side)
  //
  Double_t X1;
  Double_t Y1 = 24.3;
  Double_t Y2 = 0;
  Double_t Z1 = 0;
  for (Int_t i=0; i<8; i++)
  {
    Double_t sign = 1.;
    Double_t X2;
    if (i==1||i==2||i==5||i==6) sign = -1;
    X2 = aX_PMT[i];
    if (aY_PMT[i]>0) { Y1 =  24.3; Y2 = aY_PMT[i]-fPMTdz/2.; }
    else             { Y1 = -24.3; Y2 = aY_PMT[i]+fPMTdz/2.; }
    // if (i<4) Z1 = -kADACelldz/2. -0.1;
    if (i<4) Z1 = +kADACelldz/2. +0.1;
    else     Z1 = +kADACelldz/2. +0.1; 
    TGeoVolume * voFiber_ADA = 0;
    fX1FiberShort[i] = sign*fX_ADA_WLS_s;
    fX1FiberLong [i] = sign*fX_ADA_WLS_l;
    fX2Fiber[i] = X2;
    fY1Fiber[i] = Y1;
    fY2Fiber[i] = Y2;
    voFiber_ADA = new TGeoVolume(Form("voFiberADAShort_%d",i), MakeFiberBundle(Form("FiberADAShort_%d",i), fX1FiberShort[i], X2, Y1, Y2, Z1), medADASci);
    voFiber_ADA -> SetLineColor(kCyan);
    voADAPMTarray->AddNode (voFiber_ADA, 1); 
    voFiber_ADA = new TGeoVolume(Form("voFiberADALong_%d",i) , MakeFiberBundle(Form("FiberADALong_%d",i) , fX1FiberLong [i], X2, Y1, Y2, Z1), medADASci);
    voFiber_ADA -> SetLineColor(kCyan);
    voADAPMTarray->AddNode (voFiber_ADA, 1); 
  }
  printf( " [ADA] %14s %14s %14s %14s %14s\n", "fX1FiberShort", "fX1FiberLong", "fX2Fiber", "fY1Fiber", "fY2Fiber" );
  for (Int_t i=0; i<8; i++) {
    printf( " [%3d] %14f %14f %14f %14f %14f\n", i, fX1FiberShort[i] , fX1FiberLong [i] , fX2Fiber[i] , fY1Fiber[i] , fY2Fiber[i] );

  }

  new TGeoBBox("shADAWallPlate1", 25., 25., 0.15);
  new TGeoTube("shADAWallPlate2",  0.,  7., 0.30);
  TGeoCompositeShape * shADAWallPlate = new TGeoCompositeShape("shADAWallPlate", "shADAWallPlate1-shADAWallPlate2");
  TGeoVolume * voADAWallPlate = new TGeoVolume("voADAWallPlate", shADAWallPlate, kMedSteelSh);
  // TGeoVolume * voADAWallPlate = new TGeoVolume("voADAWallPlate", new TGeoBBox(25., 25., 0.15), kMedSteelSh);
  voADAWallPlate -> SetLineColor(kGray+3);
  ad->AddNode ( voADAWallPlate  , 1 , new TGeoTranslation ( 0. , 0. , 1701.5                      ) ); 
  ad->AddNode ( vADAarray       , 1 , new TGeoTranslation ( 0. , 0. , kPosADA - kADACelldz/2. -0.1) ); 
  ad->AddNode ( vADAarray       , 2 , new TGeoTranslation ( 0. , 0. , kPosADA + kADACelldz/2. +0.1) ); 
  ad->AddNode ( voADsupport     , 1 , new TGeoTranslation ( 0. , 0. , kPosADA                     ) ); 
  ad->AddNode ( voADsupport     , 2 , new TGeoCombiTrans  ( 0. , 0. , kPosADA,  Rz180             ) ); 
  ad->AddNode ( voADAPMTarray   , 1 , new TGeoTranslation ( 0. , 0. , kPosADA                     ) ); 
  // FAKE-VOL
  //--------------------------------------------------------------------------------//
  // ad->AddNode ( voFakeVol       , 1 , new TGeoTranslation ( 0. , 0. , 1685.0) );  // before A-Side AD
  // ad->AddNode ( voFakeVol       , 2 , new TGeoTranslation ( 0. , 0. , 1260.0) );  // after A-Side Magnet
  // ad->AddNode ( voFakeVol       , 3 , new TGeoTranslation ( 0. , 0. ,  850.0) );  // before A-Side Magnet
  // ad->AddNode ( voFakeVol       , 4 , new TGeoTranslation ( 0. , 0. ,-1885.0) );  // after  C-Side Muon concrete shielding
  // ad->AddNode ( voFakeVol       , 5 , new TGeoTranslation ( 0. , 0. ,-1948.0) );  // after  C-Side Trigger Shield
  // ad->AddNode ( voFakeVol       , 6 , new TGeoTranslation ( 0. , 0. ,-1951.5) );  // before C-Side Trigger Shield
  //--------------------------------------------------------------------------------//

  // ==========================================================================
  //
  // Define ADC (2014, May 4) Updated 2015, Jan 22
  //
  // ==========================================================================

  /////////////////////////////////////////////////////////////////////////////
  /// ADC in the tunnel                                                     ///
  /////////////////////////////////////////////////////////////////////////////
  new TGeoBBox( "shADCbox" , kADCCellSideX/2.0, kADCCellSideY/2.0, kADCCelldz/2.0 );
  new TGeoTube( "shADCHole",               0. , kADCBeamPipeR    , kADCCelldz     );
  X = kShiftX + kADCCellSideX * 0.5;
  Y = kShiftY + kADCCellSideY * 0.5;
  ( new TGeoTranslation("trADCbox", X, Y, 0.) ) -> RegisterYourself();
  // 
  TGeoVolume * vADCpad = new TGeoVolume( "ADCpad", 
    new TGeoCompositeShape("shADCpad", "shADCbox:trADCbox-shADCHole"), medADCSci );      
  vADCpad->SetLineColor( kColorADC );
  
  /// Creating Sector for Tunnel (Asembly:  Scintillator Pad + Light guide + PM )
  TGeoVolume *voADC  = new TGeoVolumeAssembly("ADCsec");
  // Add PAD
  voADC->AddNode( vADCpad, 1, 0);
  // Add ADC WLS Short bar
  voADC->AddNode( vADC_WLS_s, 1, 
      new TGeoTranslation( 0.1 + WLS_dx/2.0, kADCBeamPipeR + WLS_SideC_Short_dy/2.0, 0.0) ); 
  // Add ADC WLS Long  bar
  voADC->AddNode( vADC_WLS_l, 1, 
      new TGeoTranslation( 0.04 + WLS_dx/2.0 + kADCCellSideX + kShiftX, kShiftY + WLS_SideC_Long_dy/2.0, 0.0) ); 
  
  /// Assembling ADC adding the 4 sectors                 //  Sectors
  TGeoVolume *vADCarray = new TGeoVolumeAssembly("ADC");  //        ^ y
  vADCarray->AddNode( voADC, 1 );                         //        |   
  vADCarray->AddNode( voADC, 2, Ry180 );                  //   2    |   1
  vADCarray->AddNode( voADC, 3, Rz180 );                  // --------------->  x  
  vADCarray->AddNode( voADC, 4, Rx180 );                  //   3    |   4
                                                          //        |
                                                                             

  // ==========================================================================
  //
  // Add ADC to AD volume
  //
  // Note to future maintainers: 
  // In previous AliRoot versions the position z = -1900.75 corresponded 
  // to the end of the YSAA3_CC_BLOCK (concrete block shielding just before 
  // the C-Side LHC wall). Now this has been fixed to agree with reality. 
  // The YSAA3_CC_BLOCK starts at 1800.75 and ends at 1880.75 cm.
  //
  // Ernesto Calvo and Alberto Gago.
  // - ecalvovi@cern.ch
  // - agago@pucp.edu.pe
  //
  // ==========================================================================
  
  // *  ad -> AddNode(vADCarray , 1, new TGeoTranslation(0., 0., kZendADC2 + kADCCelldz/2.)); // Tunnel
  // *  ad -> AddNode(vADCarray , 2, new TGeoTranslation(0., 0., kZbegADC1 - kADCCelldz/2.)); // Tunnel
  // *  const Float_t kZbegADC1 = -kZbegFrontBar-2.
  // *  const Float_t kZendADC2 = -1959.0;            // (ecalvovi@cern.ch) 
  
  switch (fADCPosition ) {
    case kADCInTunnel:
      {
        // const Float_t kZbegADC1 = -kZbegFrontBar-2.;  // (ecalvovi@cern.ch) 
        // const Float_t kZendADC2 = -1959.0;            // (ecalvovi@cern.ch) 
        // ad -> AddNode(vADCarray , 1, new TGeoTranslation(0., 0., kZendADC2 + kADCCelldz/2.)); // Tunnel
        // ad -> AddNode(vADCarray , 2, new TGeoTranslation(0., 0., kZbegADC1 - kADCCelldz/2.)); // Tunnel
        const Float_t kPosADC = -kZbegFrontBar-2.-3.0-0.3;  // 3.0 = (5.6 + 0.2 + 0.2)/2. // (ecalvovi@cern.ch) 
        printf("CreateAD: kPosADC=%8.2f\n", kPosADC);
        ad -> AddNode(vADCarray,   1, new TGeoTranslation(0., 0., kPosADC - kADCCelldz/2. - 0.23)); // Tunnel // ADC1
        ad -> AddNode(vADCarray,   2, new TGeoTranslation(0., 0., kPosADC + kADCCelldz/2. + 0.23)); // Tunnel // ADC2
        ad -> AddNode(voADsupport, 3, new TGeoTranslation(0., 0., kPosADC));
        ad -> AddNode(voADsupport, 4, new TGeoCombiTrans (0., 0., kPosADC, Rz180));
        break;
      }
    case kADCInCavern:
      {
        printf("FATAL: vADCInCavern is now obsolete!");
        exit(1);
        break;
      }
    case kADCInBoth:
      {
        printf("FATAL: vADCInBoth   is now obsolete!");
        exit(1);
        break;
      }
  }


  // ==========================================================================
  // 
  // Add structure volumes to vADCstruct and vADAstruct volume assemblies
  //
  // ==========================================================================
  // 
  // Add some colors:
  //
  voWallBigPlate -> SetLineColor(kGray);
  voWallSqrPlate -> SetLineColor(kGray);
  voSaa3Rod      -> SetLineColor(kGray);
  voSaa3EndPlate -> SetLineColor(kGray);
  voFrontBar     -> SetLineColor(kGray);
  voLatBar       -> SetLineColor(kGray);
  voCoil         -> SetLineColor(kOrange+7);
  voValveVA      -> SetLineColor(kAzure+1);
  voValveVB      -> SetLineColor(kAzure+1);
  voValveVC      -> SetLineColor(kAzure+1);
  voValveVD      -> SetLineColor(kAzure+1);

  // 
  // Create Grid of U-Profiles Bars near the Wall
  //
  const Int_t nvp = 5;
  TGeoVolume * uprofiles[nvp];
  uprofiles[0] = Make_UProfile("AD_UProfileV_1", 295.0, kMedSteelSh, 10, 5, 0.85, 0.6);
  uprofiles[1] = Make_UProfile("AD_UProfileV_2", 457.5, kMedSteelSh, 10, 5, 0.85, 0.6);
  uprofiles[2] = Make_UProfile("AD_UProfileV_3", 472.5, kMedSteelSh, 10, 5, 0.85, 0.6);
  uprofiles[3] = Make_UProfile("AD_UProfileV_4", 391.0, kMedSteelSh, 10, 5, 0.85, 0.6);
  uprofiles[4] = Make_UProfile("AD_UProfileV_5", 295.0, kMedSteelSh, 10, 5, 0.85, 0.6);

  TGeoVolumeAssembly * voADC_UProfileGrid = new TGeoVolumeAssembly("voADC_UProfileGrid");
  Double_t xv[nvp] = {-119.8 , -19.8 , 80.2 , 180.2 , 280.2 }; 
  Double_t yv[nvp] = {  43.1 ,   9.4 , 34.6 ,  43.1 ,  43.1 }; 
  z = -kZwall+dAlWallThick;
  for (Int_t i=0; i<nvp; i++)
  {
    voADC_UProfileGrid->AddNode(uprofiles[i], 1, new TGeoCombiTrans(xv[i],yv[i], z, Rx90));
  }
  // 
  // Create Grid of U-Profiles Bars near the Wall
  //
  TGeoVolume * volProfH = Make_UProfileH("AD_UProfileH", kMedSteelSh);
  Double_t xh[nvp] = {-69.8, 30.2, 130.2, 230.2};
  Double_t yh[nvp] = {-64.6, 35.4, 135.4, 235.4};

  TGeoRotation * rotXY90 = NULL;
  rotXY90 = (TGeoRotation*) Rx90 -> MakeClone();
  rotXY90 -> MultiplyBy(Ry90, 1);

  for (Int_t iy=0; iy<4; iy++)
  {
    for (Int_t ix=0; ix<4; ix++)
    {
      if (iy==3) if (ix==0 || ix==3) continue;
      voADC_UProfileGrid->AddNode(volProfH, iy*4 + ix, new TGeoCombiTrans(xh[ix],yh[iy], z, rotXY90));
    }
  }
      
  TGeoVolumeAssembly * vADCstruct = new TGeoVolumeAssembly("voADCStruct");
  vADCstruct->AddNode(voSaa3EndPlate, 1, new TGeoTranslation( 0., 0., kZendAbs + 1.95/2.));
  z = kZwall;
  vADCstruct->AddNode(voWallBigPlate, 1, new TGeoTranslation(0., 0., z - 0.5 * dAlWallThick ));
  vADCstruct->AddNode(voWallSqrPlate, 1, new TGeoTranslation(0., 0., z + 0.5 * dAlWallThick ));
  z = kZendAbs + 1.95 + dzRodL/2.; 
  vADCstruct->AddNode(voSaa3Rod,  1, new TGeoTranslation(  12.5, -12.75, z));
  vADCstruct->AddNode(voSaa3Rod,  2, new TGeoTranslation(  12.5,  12.75, z));
  vADCstruct->AddNode(voSaa3Rod,  3, new TGeoTranslation( -12.5, -12.75, z));
  vADCstruct->AddNode(voSaa3Rod,  4, new TGeoTranslation( -12.5,  12.75, z));
  vADCstruct->AddNode(voValve,    1, new TGeoTranslation( 0., 0., zPosValve));
  vADCstruct->AddNode(voMoVMAOI,  1, new TGeoTranslation( 0., 0., kZbegVMAOI));
  vADCstruct->AddNode(voFrontBar, 1, new TGeoTranslation( 0., 0., kZbegFrontBar + dzF/2.));
  z = kZbegCoil;
  vADCstruct->AddNode(voCoil, 1, new TGeoCombiTrans(  3.6 + dz/2., 0., z, Ry90m));
  vADCstruct->AddNode(voCoil, 2, new TGeoCombiTrans(  3.6 + dz/2., 0., z, new TGeoRotation((*Ry90m)*(*Rx180))));
  vADCstruct->AddNode(voCoil, 3, new TGeoCombiTrans( -3.6 - dz/2., 0., z, Ry90m));
  vADCstruct->AddNode(voCoil, 4, new TGeoCombiTrans( -3.6 - dz/2., 0., z, new TGeoRotation((*Ry90m)*(*Rx180))));
  z = kZbegFrontBar + dzF + kdzLatBar/2.;
  vADCstruct->AddNode(voLatBar, 1, new TGeoTranslation(  31.9, 0., z));
  vADCstruct->AddNode(voLatBar, 2, new TGeoTranslation( -31.9, 0., z));
  vADCstruct->AddNode(CreatePmtBoxC(), 1); 
  // Create C-Side Fiber Bundles:
  // CreateCurvedBundles(ad);
  //
  // Color for voADsupp17 
  voADsupp17MainProp ->SetLineColor(kGray+1);
  voADsupp17Base     ->SetLineColor(kGray+1);
  voADsuppTopBracket ->SetLineColor(kGray+1);
  // Color for voADsupp18 
  voADsupp18MainProp -> SetLineColor(kGray+1);
  // Color for voADsuppIBeam
  voADsuppIBeam -> SetLineColor(kGreen+3);
  voADsuppIBeamV -> SetLineColor(kGreen+3);
  // voADsuppIBeam -> SetLineColorAlpha(kGreen+3, 0.);
  // voADsuppIBeam -> SetFillColorAlpha(kGreen+3, 0.5);
  // voADsuppIBeam -> SetTransparency(16);
  //
  TGeoVolume * voBLMsupport = gGeoManager->MakeBox("voBLMsupport", kMedAlu, 1.5/2.0, 30.0/2.0, 6.0/2.0);
  voBLMsupport -> SetLineColor(kGray);
  //
  // TGeoVolume         * voSupportZEM           = CreateSupportZEM(); 
  // TGeoVolume         * voPump_VAMQF           = CreatePump();
  // TGeoVolume         * voMagnet               = CreateMagnetYoke(); 
  // TGeoVolumeAssembly * voBLM                  = CreateBLM(); 
  // TGeoVolumeAssembly * voPipeOvalyzed         = CreatePipeOvalyzed(kMedSteelSh);
  TGeoVolumeAssembly * voVacuumChamberSupport = CreateVacuumChamberSupport();
  TGeoVolumeAssembly * voShield = CreateADAShielding();
  //
  TGeoVolumeAssembly * vADAstruct = new TGeoVolumeAssembly("voADAStruct");
  vADAstruct->AddNode(voADACablingVBar      , 1 , new TGeoTranslation(28.25 , 0.00  , 1604.00 )); 
  vADAstruct->AddNode(voADACablingVBar      , 2 , new TGeoTranslation(98.25 , 0.00  , 1604.00 )); 
  // FIXME: vADAstruct->AddNode(voADACablingHBar      , 1 , new TGeoTranslation(63.25 , 0.00  , 1604.00 )); 
  vADAstruct->AddNode(voADACablingHBar      , 2 , new TGeoTranslation(63.25 , 0.00  , 1604.00 )); 
  vADAstruct->AddNode(voADAMagnetCableArray , 1 , new TGeoTranslation(26.25 , 0.00  , 1608.00 )); 
  vADAstruct->AddNode(voADsupp17            , 1 , new TGeoTranslation(   0. , -36.5 , 1541.   )); 
  vADAstruct->AddNode(voADsupp18            , 1 , new TGeoTranslation(   0. , -36.5 , 1612.   )); // As Measured
  vADAstruct->AddNode(voADsupp18            , 2 , new TGeoTranslation(   0. , -36.5 , 1306.8  )); // As Measured
  // vADAstruct->AddNode(voADsupp18            , 2 , new TGeoTranslation(   0. , -36.5 , 1295.8  )); // As LHC Drawing
  vADAstruct->AddNode(voADsuppIBeam         , 1 , new TGeoTranslation(   0. , -36.5 , 1701.   ));
  vADAstruct->AddNode(voADsuppIBeamV        , 1 , new TGeoCombiTrans (   0. , -36.5-9.6-(81.8-9.6)/2.0,         1295.47+9.6/2.0, Rx90   )  );
  vADAstruct->AddNode(voADsuppIBeamV        , 2 , new TGeoCombiTrans (   0. , -36.5-9.6-(81.8-9.6)/2.0, 294.0 + 1295.47+9.6/2.0, Rx90   )  );
  vADAstruct->AddNode(voVacuumChamberSupport, 1 , new TGeoTranslation(   0. ,   0.0 , 1075.+125.         ));
  vADAstruct->AddNode(voVacuumChamberSupport, 2 , new TGeoCombiTrans (   0. ,   0.0 , 1075.-125. , Ry180 ));
  vADAstruct->AddNode(voShield              , 1 , new TGeoTranslation(   0. ,   0.0 , 1665.5             ));
  // vADAstruct->AddNode(voPipeOvalyzed        , 1 , new TGeoTranslation(   0. ,   0.0 , 1075.47 )  );
  // vADAstruct->AddNode(voBLM                 , 1 , new TGeoTranslation(  12.5,   0.0 , 1459.80 )  );
  // vADAstruct->AddNode(voBLMsupport          , 1 , new TGeoTranslation(  12.5, -21.5 , 1459.80 )  );
  // vADAstruct->AddNode(voMagnet              , 1 , new TGeoTranslation(   0. ,   0.0 , 1075.00 )  );
  // vADAstruct->AddNode(voSupportZEM          , 1 , new TGeoTranslation(   0. ,   0.0 ,  804.50 )  );

  if (fADCstruct) {
    ad->AddNode(voADC_UProfileGrid, 1       );
    ad->AddNode(vADCstruct        , 1, Ry180);
  }
  if (fADAstruct) { ad->AddNode(vADAstruct,1       ); }

  //
  // Add Everything to ALICE
  //
  // TGeoVolume *alice = gGeoManager->GetVolume("ALIC");
  alice->AddNode(ad, 1);
  

  // gGeoManager->DefaultColors();
  // // gGeoManager->CloseGeometry();
  // // gGeoManager->SetVisLevel(10);
  // // gGeoManager->SetVisOption(0);
  // // alice->Draw("ogl");
  /// C - Side
  return 0;
  TGLViewer * fViewer = (TGLViewer *) gPad->GetViewer3D();
  printf("TGLSAViewer: %p", fViewer);
  // from https://root.cern.ch/phpBB3/viewtopic.php?f=3&t=8143
  Double_t zoom      = 3.98;
  Double_t fov       = 2.1; // field of view
  Double_t dolly     =   30;  // distance from center
  Double_t center[3] = {0.,0.,-1670.};
  Double_t vRotate   =  65.;
  Double_t hRotate   = -10.;

  TGLViewer::ECameraType camera = TGLViewer::kCameraPerspXOZ;
  fViewer->SetCurrentCamera(camera);
  fViewer->CurrentCamera().SetExternalCenter(1);
  fViewer->CurrentCamera().SetCenterVecWarp(0., 0., 1075);
  fViewer->SetPerspectiveCamera (camera, fov, dolly, center, hRotate * TMath::DegToRad(), vRotate * TMath::DegToRad());
  /// A - Side
  ////////// // You can call gPad->GetViewer3D() and cast this to TGLViewer (or TGLSAViewer)
  ////////// // TGLSAViewer * fViewer = (TGLSAViewer *) gPad->GetViewer3D();
  ////////// TGLViewer * fViewer = (TGLViewer *) gPad->GetViewer3D();
  ////////// printf("TGLSAViewer: %p", fViewer);
  ////////// fViewer->CurrentCamera().SetExternalCenter(1);
  ////////// fViewer->CurrentCamera().SetCenterVecWarp(0., 0., 1690);
  ////////// // from https://root.cern.ch/phpBB3/viewtopic.php?f=3&t=8143
  ////////// Double_t zoom      = 0.98;
  ////////// Double_t fov       = 3.5; // field of view
  ////////// Double_t dolly     = 3;  // distance from center
  ////////// Double_t center[3] = {0,0,1690};
  ////////// Double_t vRotate   = 75.;
  ////////// Double_t hRotate   = 0;

  ////////// TGLViewer::ECameraType camera = TGLViewer::kCameraPerspXOZ;
  ////////// fViewer->SetCurrentCamera(camera);
  ////////// fViewer->SetPerspectiveCamera (camera, fov, dolly, center, hRotate * TMath::DegToRad(), vRotate * TMath::DegToRad());
  printf("<=== AliADv1::CreateAD(): ver=[Feb 3st, 2015]; contact=[ecalvovi@cern.ch]\n");
  return 0; 
  //
  // TGLViewer::ECameraType camera = TGLViewer::kCameraOrthoXOZ;
  // TGLViewer *v = (TGLViewer *) gPad->GetViewer3D();
  // v->SetCurrentCamera(camera);
  // v->CurrentCamera().SetExternalCenter(kTRUE);

  // Double_t zoom      = 0.98;
  // Double_t dolly     = 4;
  // Double_t center[3] = {0,0,1690};
  // Double_t hRotate   = 0.01;
  // Double_t vRotate   = 0.01;

  // v->SetOrthoCamera(camera, zoom, dolly, center, hRotate, vRotate);
  // v->DoDraw();
}


/*
TGeoVolumeAssembly * CreateAD()
{
  printf("===> AliADv1::CreateAD(): ver=[Feb 3st, 2015]; contact=[ecalvovi@cern.ch]\n");
  Int_t nvertices;
  //
  // Define Rotations used
  //
  Rx90m = new TGeoRotation("Rx90m",   0., -90.,   0.) ;
  Rx90  = new TGeoRotation("Rx90" ,   0.,  90.,   0.) ;
  Rx180 = new TGeoRotation("Rx180",   0., 180.,   0.) ;  //   4    |   1
  Rz180 = new TGeoRotation("Rz180", 180.,   0.,   0.) ;  // --------------->  x  
  Ry180 = new TGeoRotation("Ry180", 180., 180.,   0.) ;  //   3    |   2
  Ry90m = new TGeoRotation("Ry90m",  90., -90., -90.) ;
  Ry90  = new TGeoRotation("Ry90" ,  90.,  90., -90.) ;
  // Get Mediums needed. 
  // TGeoMedium * kMedAlu       = gGeoManager->GetMedium("AD_Alum");   // Aluminium 
  // TGeoMedium * kMedSteelSh   = gGeoManager->GetMedium("AD_ST_C0");  // Stainless Steel 
  // TGeoMedium * kMedVacuum    = gGeoManager->GetMedium("AD_VA_C0");  // Stainless Steel 
  // TGeoMedium * kMedCopper    = gGeoManager->GetMedium("AD_Cu_C0");  // Stainless Steel 
  // TGeoMedium * kMedPVC       = gGeoManager->GetMedium("AD_PVC"  );  // Stainless Steel 
  /////////////////////////////////////////////////////////////////////////////
  //
  // Add Cables for A-Side
  // No LHC drawing found!
  //
  new TGeoBBox("shADACablingVBar1", 2.0, 100.0, 2.0);
  new TGeoBBox("shADACablingVBar0", 1.8, 102.0, 1.8);
  TGeoVolume* voADACablingVBar = new TGeoVolume("voADACablingVBar", 
    new TGeoCompositeShape("shADACablingVBar", "shADACablingVBar1-shADACablingVBar0"), kMedSteelSh);
  new TGeoBBox("shADACablingHBar1",  70.0, 2.0, 2.0);
  new TGeoBBox("shADACablingHBar0",  72.0, 1.8, 1.8);
  TGeoVolume* voADACablingHBar = new TGeoVolume("voADACablingHBar", 
    new TGeoCompositeShape("shADACablingHBar", "shADACablingHBar1-shADACablingHBar0"), kMedSteelSh);
  voADACablingVBar->SetLineColor(kGray+1);
  voADACablingHBar->SetLineColor(kGray+1);

  /////////////////////////////////////////////////////////////////////////////
  // 
  // Vertical Cables A-Side
  // No LHC drawing found!
  // There are ~ 37 cables in a dx space of 74 cm. That makes cables of 2 cm diameter
  // No info about the cable. Taking as an aproximation a cable from Riyadh catalog:
  // model: 000101xx15 
  // nominal cross section: 1 x 120 
  // number of wires in conductor: 37 
  // diameter of conductor: 14.21 
  // Insulation thickness: 1.6 
  // Overall diameter: 18.8
  //
  // It has 37 17.64
  //
  Double_t rCable     = 2.0/2.;
  Double_t rCableCore = 1.4/2.; 
  Double_t lCable     = 200.;
  // TGeoTube   * shADCableExt  = new TGeoTube   ("shADCableExt" , rCableCore, rCable, lCable/2.);
  // TGeoTube   * shADCableCore = new TGeoTube   ("shADCableCore",         0., rCable, lCable/2.);
  TGeoVolume * voADCableExt  = new TGeoVolume ( "voADCableExt"  , new TGeoTube ( "shADCableExt"  , rCableCore+0.01 , rCable          , lCable/2. )  , kMedPVC) ;
  TGeoVolume * voADCableCore = new TGeoVolume ( "voADCableCore" , new TGeoTube ( "shADCableCore" , 0.              , rCableCore-0.01 , lCable/2. )  , kMedCopper);
  voADCableExt  -> SetLineColor(kGray+3);
  voADCableCore -> SetLineColor(kOrange+2);
  TGeoVolumeAssembly * voADAMagnetCable = new TGeoVolumeAssembly("voADAMagnetCable");
  voADAMagnetCable->AddNode(voADCableExt , 1, Rx90);
  voADAMagnetCable->AddNode(voADCableCore, 1, Rx90);




  TGeoVolumeAssembly * voADAMagnetCableArray= new TGeoVolumeAssembly("voADAMagnetCableArray");
  
  for (Int_t i=0; i<37; i++)
  {
    voADAMagnetCableArray -> AddNode(voADAMagnetCable, i+1, new TGeoTranslation(+rCable + 2.*rCable*i,0,0));
  }
  // Make I-Beam support structure in RB24/3 (A-Side)
  // LHC Drawing: LHCVC2U_0034
  // Dimensions are approximated.
  // Probably a type "HE 100 A" Wide Flange Beam.
  // flange web      :  0.5 cm
  // flange width    : 10.0 cm
  // flange height   :  9.6 cm
  // flange thickness:  0.8 cm
  TGeoXtru * shADsuppIBeam = new TGeoXtru(2);
  shADsuppIBeam->SetNameTitle("shADsuppIBeam","shADsuppIBeam");
  Double_t HIBeam, hIBeam, h2, W, w;
  HIBeam  = 9.6/2.;
  hIBeam  = 0.8;
  h2 = hIBeam;
  w  =  0.5/2.;
  W  = 10.5/2.;
  Double_t xIBeam[] = { -W      , -W            , -w        , -w        , -W            , -W     , W      , W             , w         , w         , W             , W       }; 
  Double_t yIBeam[] = { -HIBeam , hIBeam-HIBeam , h2-HIBeam , HIBeam-h2 , HIBeam-hIBeam , HIBeam , HIBeam , HIBeam-hIBeam , HIBeam-h2 , h2-HIBeam , hIBeam-HIBeam , -HIBeam }; 
  nvertices = sizeof(xIBeam)/sizeof(Double_t);
  shADsuppIBeam->DefinePolygon(nvertices,xIBeam,yIBeam);
  shADsuppIBeam->DefineSection(0, -430., 0., -HIBeam, 1.0); // index, Z position, offset (x,y) and scale for first section
  shADsuppIBeam->DefineSection(1,    0., 0., -HIBeam, 1.0); // idem, second section
  Double_t fOriginBeamCut1[3] = { 3.25, -HIBeam, -4.5};
  Double_t fOriginBeamCut2[3] = {-3.25, -HIBeam, -4.5};
  new TGeoBBox("shADBeamCut1", W/2., HIBeam + hIBeam, 3., fOriginBeamCut1);
  new TGeoBBox("shADBeamCut2", W/2., HIBeam + hIBeam, 3., fOriginBeamCut2);
  TGeoCompositeShape * shADsuppIBeamCutted = new TGeoCompositeShape("shADsuppIBeamCutted", "(shADsuppIBeam-shADBeamCut1)-shADBeamCut2");
  TGeoVolume * voADsuppIBeam = new TGeoVolume("voADsuppIBeam", shADsuppIBeamCutted, kMedSteelSh);

  // Vertical I-Beam:
  TGeoVolume * voADsuppIBeamV = MakeVolIBeam("voADsuppIBeamV"    , kMedSteelSh, 10.5, 9.6, 0.5, 1.6,81.8-9.6);

  /////////////////////////////////////////////////////////////////////////////
  //
  // Make Top Bracket, part of 
  // Sliding support (structure 18 in RB24/3)
  // LHC Drawing: LHCVC2U_0026 (item 5)
  // 
  new TGeoBBox("shADsupp1Box",                  11./2., 15./2., 1. );
  new TGeoBBox("shADsupp1HoleSqr",               8./2.,  4./2., 1.2);  // larger y (4 instead of 3) just in case
  new TGeoTube("shADsupp1HoleCircle",               0.,    0.5, 1.2);
  new TGeoBBox("shADsupp1HoleRectangle",            1.,    0.5, 1.2);
  TGeoArb8 * shADsupp1HoleTrg = new TGeoArb8("shADsupp1HoleTrg", 1.2);
  Float_t trgX[] = { -2.5, 0.0, 2.5,  0.0 };
  Float_t trgY[] = {  0.0, 2.5, 0.0, -2.5 };
  for (Int_t i=0; i<4; i++) 
  {
    shADsupp1HoleTrg -> SetVertex(i,   trgX[i], trgY[i]);
    shADsupp1HoleTrg -> SetVertex(i+4, trgX[i], trgY[i]);
  }
  (new TGeoTranslation("trADsupp1Box",              0., 15./2., 0. ))->RegisterYourself();
  (new TGeoTranslation("trADsupp1HoleSqr",          0.,    14., 0. ))->RegisterYourself();
  (new TGeoTranslation("trADsupp1HoleRndRecL1",     3.,    2.5, 0. ))->RegisterYourself();
  (new TGeoTranslation("trADsupp1HoleRndRecL2",     3.,    5.5, 0. ))->RegisterYourself();
  (new TGeoTranslation("trADsupp1HoleRndRecL3",     3.,    8.5, 0. ))->RegisterYourself();
  (new TGeoTranslation("trADsupp1HoleRndRecR1",    -3.,    2.5, 0. ))->RegisterYourself();
  (new TGeoTranslation("trADsupp1HoleRndRecR2",    -3.,    5.5, 0. ))->RegisterYourself();
  (new TGeoTranslation("trADsupp1HoleRndRecR3",    -3.,    8.5, 0. ))->RegisterYourself();
  (new TGeoTranslation("trADsupp1HoleCircleL",      1.,     0., 0. ))->RegisterYourself();
  (new TGeoTranslation("trADsupp1HoleCircleR",     -1.,     0., 0. ))->RegisterYourself();
  new TGeoCompositeShape("shADsupp1HoleRoundedRectangle", 
      "shADsupp1HoleCircle:trADsupp1HoleCircleL + "
      "shADsupp1HoleCircle:trADsupp1HoleCircleR + "
      "shADsupp1HoleRectangle");
  TGeoCompositeShape * shADsuppTopBracket = new TGeoCompositeShape("shADsuppTopBracket", 
      "shADsupp1Box:trADsupp1Box - "
      "shADsupp1HoleRoundedRectangle:trADsupp1HoleRndRecL1 - "
      "shADsupp1HoleRoundedRectangle:trADsupp1HoleRndRecL2 - "
      "shADsupp1HoleRoundedRectangle:trADsupp1HoleRndRecL3 - "
      "shADsupp1HoleRoundedRectangle:trADsupp1HoleRndRecR1 - "
      "shADsupp1HoleRoundedRectangle:trADsupp1HoleRndRecR2 - "
      "shADsupp1HoleRoundedRectangle:trADsupp1HoleRndRecR3 - "
      "shADsupp1HoleSqr:trADsupp1HoleSqr - "
      "shADsupp1HoleTrg");
  TGeoVolume * voADsuppTopBracket = new TGeoVolume("voADsuppTopBracket", shADsuppTopBracket, kMedAlu);
  

  /////////////////////////////////////////////////////////////////////////////
  //
  // Make support BASE for structure 17 and 18 in RB24/3 (A-Side)
  // LHC Drawing: LHCVC2U_0024
  
  new TGeoBBox("shADsupp17BaseBox",   16./2., 0.65/2., 33./2. );
  new TGeoBBox("shADsupp17VertXBox",  16./2.,   8./2., 0.65/2. );
  new TGeoBBox("shADsupp17VertZBox",  1.3/2.,   8./2.,  25./2. );
  (new TGeoTranslation("trADsupp17BaseBox",  0., 0.65/2., -12.15 ))->RegisterYourself();
  (new TGeoTranslation("trADsupp17VertXBox", 0.,  8.0/2., -12.15 +  8.5 + 0.65/2. ))->RegisterYourself();
  (new TGeoTranslation("trADsupp17VertZBox", 0.,  8.0/2., -12.15 -4. ))->RegisterYourself();
  TGeoCompositeShape * shADsupp17Base = new TGeoCompositeShape("shADsupp17Base", 
      "shADsupp17BaseBox:trADsupp17BaseBox + "
      "shADsupp17VertXBox:trADsupp17VertXBox + "
      "shADsupp17VertZBox:trADsupp17VertZBox");
  TGeoVolume * voADsupp17Base = new TGeoVolume("voADsupp17Base", shADsupp17Base, kMedAlu);
  /////////////////////////////////////////////////////////////////////////////
  //
  // Make support structure 17 in RB24/3 (A-Side)
  // Fixed support
  // LHC Drawing: LHCVC2U_0028
  // LHC Drawing: LHCVC2U_0032
  
  // Main Prop
  TGeoXtru * shADsupp17MainPropLat = new TGeoXtru(2);
  shADsupp17MainPropLat->SetNameTitle("shADsupp17MainPropLat","shADsupp17MainPropLat");
  Double_t x17[] = {  0.0 , 3 , 3.  , 3.   , 3.6  , 5.   , 5.  , 0. };
  Double_t y17[] = {  0.0 , 0 , 0.8 , 7.55 , 7.55 , 8.95 , 28. , 28.};
  nvertices = sizeof(x17)/sizeof(Double_t);
  shADsupp17MainPropLat->DefinePolygon(nvertices,x17,y17);
  shADsupp17MainPropLat->DefineSection(0, 0.00, 0., 0., 1.0); // Z position, offset and scale for first section
  shADsupp17MainPropLat->DefineSection(1, 0.65, 0., 0., 1.0); // idem, second section
  (new TGeoCombiTrans("ctADsupp17BaseLatL"      , 5. - 0.65 , 0.65          , 0. , Ry90 ))->RegisterYourself();
  (new TGeoCombiTrans("ctADsupp17BaseLatR"      , -5.       , 0.65          , 0. , Ry90 ))->RegisterYourself();
  (new TGeoTranslation("trADsupp17MainPropVBox" , 0.        , 28./2. + 0.65 , -0.65/2.  ))->RegisterYourself();
  new TGeoBBox("shADsupp17MainPropVBox", 10./2., 28./2., 0.65/2.);
  TGeoCompositeShape * shADsupp17MainProp = new TGeoCompositeShape("shADsupp17MainProp",
      "shADsupp17MainPropLat:ctADsupp17BaseLatL + "
      "shADsupp17MainPropLat:ctADsupp17BaseLatR + "
      "shADsupp17MainPropVBox:trADsupp17MainPropVBox");
  TGeoVolume * voADsupp17MainProp = new TGeoVolume("voADsupp17MainProp", shADsupp17MainProp, kMedAlu);
  TGeoVolumeAssembly * voADsupp17 = new TGeoVolumeAssembly("voADsupp17");
  voADsupp17 -> AddNode(voADsupp17MainProp, 1);
  voADsupp17 -> AddNode(voADsupp17Base,     1);
  voADsupp17 -> AddNode(voADsuppTopBracket, 1, new TGeoTranslation(0., 15.90, 0.1 + 2./2.));
  



  /////////////////////////////////////////////////////////////////////////////
  //
  // Make support structure 2: part of 
  // Sliding support
  // LHC Drawing: LHCVC2U_0026
  // 
  Float_t htrp   = 4.3;  // real value??
  Float_t trpX[] = {0.  , 5.0 , 4.3,0.  };
  Float_t trpY[] = {htrp, htrp, 0. ,0.  };
  TGeoArb8 * shADsupp2trp = new TGeoArb8("shADsupp2trp", 0.6/2.);
  for (Int_t i=0; i<4; i++) 
  {
    shADsupp2trp -> SetVertex(i,   trpX[i], trpY[i]);
    shADsupp2trp -> SetVertex(i+4, trpX[i], trpY[i]);
  }
  new TGeoBBox("shADsupp2latbox", 5.0/2., (19.-htrp)/2.,  0.6/2.);
  new TGeoBBox("shADsupp2box",    0.6/2.,        19./2., 10.0/2.);
  (new TGeoTranslation("trADsupp2B"   , 0.6/2.,              19./2.,              0.))->RegisterYourself();
  (new TGeoTranslation("trADsupp2L"   , 5.0/2., (19-htrp)/2. + htrp, -10.0/2.+0.6/2.))->RegisterYourself();
  (new TGeoTranslation("trADsupp2R"   , 5.0/2., (19-htrp)/2. + htrp,  10.0/2.-0.6/2.))->RegisterYourself();
  (new TGeoTranslation("trADsupp2trpL",     0.,                  0., -10.0/2.+0.6/2.))->RegisterYourself();
  (new TGeoTranslation("trADsupp2trpR",     0.,                  0.,  10.0/2.-0.6/2.))->RegisterYourself();
  TGeoCompositeShape * shADsupp18MainProp = new TGeoCompositeShape("shADsupp18MainProp", 
      "shADsupp2trp:trADsupp2trpL + "
      "shADsupp2trp:trADsupp2trpR + "
      "shADsupp2latbox:trADsupp2L + "
      "shADsupp2latbox:trADsupp2R + "
      "shADsupp2box:trADsupp2B");
  TGeoVolume * voADsupp18MainProp = new TGeoVolume("voADsupp18MainProp", shADsupp18MainProp, kMedAlu);
  TGeoVolumeAssembly * voADsupp18 = new TGeoVolumeAssembly("voADsupp18");
  voADsupp18 -> AddNode(voADsupp18MainProp, 1, new TGeoCombiTrans (0., 9.65, 0., Ry90m));
  voADsupp18 -> AddNode(voADsupp17Base,     1, Ry180);
  voADsupp18 -> AddNode(voADsuppTopBracket, 1, new TGeoTranslation(0., 16.25, -0.1 - 2./2.));
  /////////////////////////////////////////////////////////////////////////////
  // Define Ion Pump ??
  // Drawing LHCVBU__0052
  // Vacuum - Bellows - U type
  // BODY 1 PORTS
  //
  (new TGeoCombiTrans("ctPumpVB2", 0., -13./2., 6.8-11.5/2., Rx90))->RegisterYourself();
  new TGeoTube("shIonPumpVB1o",  0.0, 10.3 /2., 11.5/2.   );
  new TGeoTube("shIonPumpVB2o",  0.0,  7.0 /2., 13.0/2.   );
  new TGeoTube("shIonPumpVB1i",  0.0, 10.0 /2., 11.5/2.+2.);
  new TGeoTube("shIonPumpVB2i",  0.0,  6.7 /2., 13.0/2.+2.);
  new TGeoCompositeShape("shIonPumpVBo", "shIonPumpVB1o+shIonPumpVB2o:ctPumpVB2");
  //
  // Continue definition of LHCVBU__0052
  //
  new TGeoCompositeShape("shIonPumpVBi", "shIonPumpVB1i+shIonPumpVB2i:ctPumpVB2");
  TGeoShape * sh3 = new TGeoCompositeShape("shIonPumpVB",  "shIonPumpVBo-shIonPumpVBi");
  TGeoVolume * voIonPumpVB = new TGeoVolume("voIonPumpVB", sh3, kMedAlu);
  // Variables 
  Double_t alpha, beta, tga2, tga, sa, ca, ctgb, d, Ro, Ri, phi1, dphi, H, L, z;
  // Drawing: LHCVSR__0054
  // Vacuum Screen - RF
  // transition flange
  alpha = 15. * TMath::DegToRad();
  beta  = 15. * TMath::DegToRad();
  tga2  = TMath::Tan(alpha*0.5);
  tga   = TMath::Tan(alpha);
  ctgb  = 1./TMath::Tan(beta );
  ca    = TMath::Cos(alpha    );
  d     = 0.3; // thickness of transition flange
  Double_t h = d/ca; // vertical distance between parallel surfaces tilted alpha degrees
  TGeoPcon * shVSRflange = new TGeoPcon("shVSRflange", 0.0, 360.0, 7);
  Ri = 9.71/2.; Ro = 11.16/2.;
  shVSRflange->DefineSection(0,         0.  ,   Ri, Ro);
  Ri = Ri - 0.25 * tga;
  shVSRflange->DefineSection(1,         0.25,   Ri, Ro);
  Ro = Ri + d*tga;
  shVSRflange->DefineSection(2,         0.25,   Ri, Ro);
  Ri = 6.3/2.; Ro = Ri + h;
  z = (9.71/2. - Ri) / tga;
  shVSRflange->DefineSection(3,            z,   Ri, Ro);
  // 
  Double_t   Dtga = 6.6*tga - 0.5*(9.71-6.3)    ;
  Double_t x = (h - 0.11 - Dtga) / (ctgb - tga) ;
  Double_t y = x * ctgb;
  z  = 6.6 - x;
  Ro = Ri + 0.11 + y;
  shVSRflange->DefineSection(4,            z,   Ri, Ro);
  z  = 6.6; 
  Ro = Ri + 0.11;
  shVSRflange->DefineSection(5,            z,   Ri, Ro);
  z  = 7.1; 
  shVSRflange->DefineSection(6,            z,   Ri, Ro);
  TGeoVolume * voVSRflange = new TGeoVolume("voVSRflange", shVSRflange, kMedAlu);
  //
  // Drawing: LHCVSR__0053
  // Vacuum Screen - RF
  // transition tube
  alpha = 10. * TMath::DegToRad();
  tga2  = TMath::Tan(alpha*0.5);
  tga   = TMath::Tan(alpha);
  ca    = TMath::Cos(alpha    );
  d     = 0.3; // thickness of Vacuum Screen RF
  Ro    = 6.7/2.; //
  Ri    = 0.;
  phi1  = 90. - 15.;
  dphi  = 30.;
  TGeoPcon * shVSR0 = new TGeoPcon("shVSR0", phi1, dphi, 6);
  shVSR0->DefineSection(0,         0.  ,   Ro-0.09, Ro);
  shVSR0->DefineSection(1,         0.45,   Ro-0.09, Ro);
  Ri=Ro-d;
  shVSR0->DefineSection(2,         0.45,      Ri, Ro);
  shVSR0->DefineSection(3, 13.37-d*tga2,      Ri, Ro);
  Ro  = Ri + d/ca;
  shVSR0->DefineSection(4,        13.37,      Ri, Ro);
  Ri += 0.63*tga;
  Ro  = Ri + d/ca;
  shVSR0->DefineSection(5,        14.0 ,      Ri, Ro);
  // printf("  Ro: %8.2f\n", Ro);
  // Make holes 
  new TGeoBBox("shHoleBody"    , 0.15, 0.60, 0.3);
  new TGeoTube("shHoleEnd", 0.  , 0.15, 0.3);
  (new TGeoTranslation("trHoleBody", 0., -0.6, 0.))->RegisterYourself();
  (new TGeoTranslation("trHoleEnd" , 0., -1.2, 0.))->RegisterYourself();
  new TGeoCompositeShape("shHole","shHoleEnd + shHoleEnd:trHoleEnd + shHoleBody:trHoleBody");
  // Single hole made. Now define some combitrans to position holes
  z = 1.3; Ro = (6.7 - d)*0.5;
  (new TGeoCombiTrans("ctHole1", 0., Ro , z, Rx90m))->RegisterYourself();
  z = 1.3 + 1*2.5;
  (new TGeoCombiTrans("ctHole2", 0., Ro , z, Rx90m))->RegisterYourself();
  z = 1.3 + 2*2.5;
  (new TGeoCombiTrans("ctHole3", 0., Ro , z, Rx90m))->RegisterYourself();
  z = 1.3 + 3*2.5;
  (new TGeoCombiTrans("ctHole4", 0., Ro , z, Rx90m))->RegisterYourself();
  // Now make a sector of RF transition tube
  new TGeoCompositeShape("shVSRsec",
   "shVSR0 - (shHole:ctHole1 + shHole:ctHole2 + shHole:ctHole3 + shHole:ctHole4)");
  // Now define rotations for each sector
  TString strSh = "shVSRsec ";
  for (Int_t i=1; i<=11; i++) {
   (new TGeoRotation(Form("rSec%d",i), 30. * i, 0. , 0.))->RegisterYourself();
   strSh+=Form("+ shVSRsec:rSec%d",i);
  }
  // printf("%s\n", strSh.Data());
  TGeoCompositeShape * shVSR = new TGeoCompositeShape("shVSR", strSh.Data());
  // Now assembly the sector to form VSR RF transition tube !
  TGeoVolume * voVSR = new TGeoVolume("voVSR", shVSR, kMedAlu);
  // 
  // Drawing: LHCVSR__0057
  // RF CONTACT
  // 
  Ro = 0.5 * 6.3;
  d  = 0.03;
  Ri = Ro - d;
  // alpha = TMath::ArcSin((7.35-1.75)/Ri);  <-- No!
  H = 0.5 * 6.9 - Ri;
  L = 28. - 7.1 + 0.45 -14. - 1.75;//7.35 - 1.75;
  // Double_t Delta = TMath::Sqrt( L*L + 4.*(H-d)*H );
  Double_t R = TMath::Sqrt((H-d)*(H-d) + L*L);
  alpha = TMath::ASin(d/R) + TMath::ASin((H-d)/R);
  //printf("alpha: %8.2f \n", alpha * TMath::RadToDeg());
  sa = TMath::Sin(alpha);
  ca = TMath::Cos(alpha);
  x = d*sa;
  y = d*ca;
  Double_t R0 =  1.75;
  Double_t R1 = 10.48;
  Double_t R2 =  0.81 + 0.28;
  phi1 = 0.; dphi = 360.; 

  TGeoPcon * shVSRcontact = new TGeoPcon("shVSRcont", phi1, dphi, 6);
  z = 0.;
  shVSRcontact->DefineSection(0, -z, Ri, Ro);
  z = R0;
  shVSRcontact->DefineSection(1, -z, Ri, Ro);
  z  += x;
  Ri += y;
  Ro  = Ri + d;
  shVSRcontact->DefineSection(2, -z, Ri, Ro);
  z  += (R1 - R0) * ca;
  Ri += (R1 - R0) * sa;
  Ro  = Ri + d;
  shVSRcontact->DefineSection(3, -z, Ri, Ro);
  // Last sections (R2)
  Double_t ab = alpha + 21. * TMath::DegToRad();
  Double_t sab = TMath::Sin(ab);
  Double_t cab = TMath::Cos(ab);
  x   = d * sab;
  y   = d/ca - d*cab;
  z  += x;
  Ri += y;
  Ro  = Ri + d/cab;
  shVSRcontact->DefineSection(4, -z, Ri, Ro);
  z  += R2 * cab;
  Ri += R2 * sab;
  Ro  = Ri + d/cab;
  shVSRcontact->DefineSection(5, -z, Ri, Ro);
  TGeoVolume * voVSRcontD = new TGeoVolume("voVSRcontD", shVSRcontact, kMedAlu);
  // Drawing: LHCVSR__0017
  // Vacuum Screen - RF
  // RF Contact flange
  phi1 = 0.;
  dphi = 360.;

  TGeoPcon * shVSRcontFlange = new TGeoPcon("shVSRcontFlange", phi1, dphi, 11);
  Ri = 6.30 /2.; Ro = 11.16 /2.; z = 0;
  shVSRcontFlange->DefineSection( 0, -z, Ri, Ro);
  Ri = 6.30 /2.; Ro = 11.16 /2.; z = 0.1;
  shVSRcontFlange->DefineSection( 1, -z, Ri, Ro);
  Ri = 6.36 /2.; Ro = 11.16 /2.; z = 0.1;
  shVSRcontFlange->DefineSection( 2, -z, Ri, Ro);
  Ri = 6.36 /2.; Ro = 11.16 /2.; z = 0.25;
  shVSRcontFlange->DefineSection( 3, -z, Ri, Ro);
  Ri = 6.36 /2.; Ro =  9.30 /2.; z = 0.25;
  shVSRcontFlange->DefineSection( 4, -z, Ri, Ro);
  Ri = 6.36 /2.; Ro =  9.30 /2.; z = 0.30;
  shVSRcontFlange->DefineSection( 5, -z, Ri, Ro);
  Ri = 6.36 /2.; Ro =  6.90 /2.; z = 0.30;
  shVSRcontFlange->DefineSection( 6, -z, Ri, Ro);
  Ri = 6.36 /2.; Ro =  6.90 /2.; z = 1.10;
  shVSRcontFlange->DefineSection( 7, -z, Ri, Ro);
  Ri = 6.36 /2.; Ro =  6.85 /2.; z = 1.15;
  shVSRcontFlange->DefineSection( 8, -z, Ri, Ro);
  Ri = 6.36 /2.; Ro =  6.50 /2.; z = 1.15;
  shVSRcontFlange->DefineSection( 9, -z, Ri, Ro);
  Ri = 6.36 /2.; Ro =  6.50 /2.; z = 1.30;
  shVSRcontFlange->DefineSection(10, -z, Ri, Ro);
  TGeoVolume * voVSRcontF = new TGeoVolume("voVSRcontF", shVSRcontFlange, kMedAlu);
  // Drawing: LHCVBU__0002
  // Bellows + End Parts
  // Vacuum - Bellows - U type

  // First make end part
  phi1 = 0. ; dphi = 360. ;
  TGeoPcon * shVBUend = new TGeoPcon("shVBUend", phi1, dphi, 6);
  Ri = 5.176; Ro = 5.4; z = 0;
  shVBUend->DefineSection( 0, z, Ri, Ro);
  Double_t dz  = 0.03;
  Ri -= dz * TMath::Tan(45. * TMath::DegToRad());
  z  += dz;
  shVBUend->DefineSection( 1, z, Ri, Ro);
  dz  = 0.04;
  Ro -= dz * TMath::Tan(15. * TMath::DegToRad());
  Ri -= dz * TMath::Tan(45. * TMath::DegToRad());
  z  += dz;
  shVBUend->DefineSection( 2, z, Ri, Ro);
  Ro  = 5.250 ;
  Ri  = 5.073;
  z   = 0.103;
  shVBUend->DefineSection( 3, z, Ri, Ro);
  Ro  = 5.25;
  Ri  = 5.02;
  z   = 0.18;
  shVBUend->DefineSection( 4, z, Ri, Ro);
  Ro  = 5.15;
  Ri  = 5.00;
  z   = 0.28;
  shVBUend->DefineSection( 5, z, Ri, Ro);

  TGeoPcon * shVBUtube26mm = new TGeoPcon("shVBUtube26mm", 0., 360., 2);
  shVBUtube26mm->DefineSection( 0, 0.28, 5.0, 5.15);
  shVBUtube26mm->DefineSection( 1, 2.60, 5.0, 5.15);
  TGeoPcon * shVBUtube9mm  = new TGeoPcon("shVBUtube9mm" , 0., 360., 2);
  shVBUtube9mm ->DefineSection( 0, 0.00, 5.0, 5.15);
  shVBUtube9mm ->DefineSection( 1, 0.62, 5.0, 5.15);
  // central part of bellow (TODO: Add plies)
  TGeoPcon * shVBUcent  = new TGeoPcon("shVBUcent" , 0., 360., 6);
  const Int_t nsec = 6;
  Double_t az [nsec] = {0.9, 0.915, 0.915, 11.885, 11.885, 11.9};
  Double_t aRi[nsec] = {5.0, 5.0  , 5.685,  5.685,  5.0  ,  5.0};
  for (Int_t i=0; i<nsec; i++) {
   z=az[i]; Ri = aRi[i]; Ro = 5.7;
   //printf("  i: %2d  z: %8.2f  Ri: %8.2f  Ro: %8.2f\n", i, z, Ri, Ro );
   shVBUcent ->DefineSection( i, z, Ri, Ro);
  }

  ( new TGeoCombiTrans("ctEnd9mm", 0., 0., 0.9, Ry180)) -> RegisterYourself();

  TGeoCompositeShape * shVBU9mm  = new TGeoCompositeShape("shVBU9mm" , "shVBUend:ctEnd9mm + shVBUtube9mm");
  TGeoCompositeShape * shVBU26mm = new TGeoCompositeShape("shVBU26mm", "shVBUend + shVBUtube26mm");

  const Int_t nsec2 = 10;
  Double_t az2 [nsec2] = {0.  , 0.43, 0.43, 0.58, 0.58 , 0.73 , 0.73, 1.05, 1.05, 1.31} ;
  Double_t aRi2[nsec2] = {5.15, 5.15, 5.03, 5.03, 5.455, 5.455, 5.03, 5.03, 5.59, 5.59} ;
  TGeoPcon * shVBUrotFlg  = new TGeoPcon("shVBUrotFlg" , 0., 360., nsec2);
  for (Int_t i=0; i<nsec2; i++) {
   z=az2[i]; Ri = aRi2[i]; Ro = 6.02;
   //printf("  i: %2d  z: %8.2f  Ri: %8.2f  Ro: %8.2f\n", i, z, Ri, Ro );
   shVBUrotFlg ->DefineSection( i, z, Ri, Ro);
  }
  TGeoVolume * voVBUrotFlg = new TGeoVolume("voVBUrotFlg", shVBUrotFlg, kMedAlu);

  // Flange 
  TGeoPcon * shVBUflg  = new TGeoPcon("shVBUflg" , 0., 360., 4);
  z  = 0; 
  Ri = 6.02; Ro = 7.6; 
  shVBUflg->DefineSection(0, z, Ri, Ro);
  z  = 1.31; 
  shVBUflg->DefineSection(1, z, Ri, Ro);
  z  = 1.31; 
  Ri = 5.15; 
  shVBUflg->DefineSection(2, z, Ri, Ro);
  z  = 2.0;
  shVBUflg->DefineSection(3, z, Ri, Ro);
  TGeoVolume * voVBUflg = new TGeoVolume("voVBUflg", shVBUflg, kMedAlu);

   
  // Add the metal plate at the end of Absorber
  // Plate:  80 cm x 80 cm x 1.95 cm (Thickness is aproximatted) (ernesto.calvo@pucp.edu.pe)
  // The End of the concrete abosorber is at kZendAbs = 1880.75

  new TGeoBBox("shBasePlate", 80./2., 80./2.,  1.95/2.);
  new TGeoTube("shHolePlate",  0.   , 12.3  ,  1.95   );
  TGeoVolume* voSaa3EndPlate  =  new TGeoVolume("voYSaa3EndPlate",
      new TGeoCompositeShape("shYSaa3EndPlate","shBasePlate-shHolePlate"),
      kMedSteelSh);
  //
  // Add Rods
  //
  // dimensions of rods
  const Double_t dzRodL = 27.0;
  const Double_t dzRodA =  4.3;
  const Double_t dzRodB =  1.3;
  new TGeoTube("shLargeRod", 0.,   1.6/2.,  dzRodL/2.);
  new TGeoTube("shRodA",     0.,   3.0/2.,  dzRodA/2.);
  new TGeoTube("shRodB",     0.,   2.3/2.,  dzRodB/2.);
  //
  ( new TGeoTranslation("trRod1", 0., 0., -dzRodL/2. + dzRodA/2.)          )->RegisterYourself();
  ( new TGeoTranslation("trRod2", 0., 0., -dzRodL/2. + dzRodA + dzRodB/2.) )->RegisterYourself();
  ( new TGeoTranslation("trRod3", 0., 0.,  dzRodL/2. - dzRodB/2.)          )->RegisterYourself();
  TGeoVolume * voSaa3Rod = new TGeoVolume("YSAA3Rod",
      new TGeoCompositeShape("shLargeRod+shRodA:trRod1 + shRodB:trRod2 + shRodB:trRod3"),
      kMedSteelSh);
  //
  // Define Valve support (VS)  (ernesto.calvo@pucp.edu.pe) 
  //
  const Double_t dyVS =  5.5; 
  const Double_t dxVS = 30.0; 
  const Double_t dzVS =  1.0;
  TGeoVolume * voVS = new TGeoVolume("voVS",
      new TGeoBBox("shVS", dxVS/2., dyVS/2., dzVS/2.),
      kMedSteelSh);
  voVS->SetLineColor(kGray);

  // Add Valve (Valve is divided in parts VA,VB,VC and VD)  (ernesto.calvo@pucp.edu.pe) 
  TGeoVolumeAssembly * voValve = new TGeoVolumeAssembly("voValve");
  //
  // Define volume VA  (ernesto.calvo@pucp.edu.pe) 
  const Double_t dxVA  = 20.3; 
  const Double_t dyVA  = 48.0; 
  const Double_t dzVA  =  6.0; // Width 
  const Double_t dz2VA =  8.5; // Full width including protuding boxes
  // Valve position  (ernesto.calvo@pucp.edu.pe) 
  const Double_t zPosValve =  kZbegValve + dz2VA/2.;
  //
  new TGeoBBox("shVAbox",       dxVA/2., dyVA/2.,     dzVA/2.);
  new TGeoBBox("shVAHbox",  -1.+dxVA/2.,   3./2.,    dz2VA/2.);
  new TGeoTube("shVAC",              0.,     7.9,    dz2VA/2.);
  new TGeoTube("shVACh",             0.,     5.3,    dz2VA   );
  // translation for shVAHbox (ernesto.calvo@pucp.edu.pe) 
  ( new TGeoTranslation("trVAH1", 0.,  12.75, 0.) )->RegisterYourself();
  ( new TGeoTranslation("trVAH2", 0., -12.75, 0.) )->RegisterYourself();

  TGeoVolume * voValveVA = new TGeoVolume("voValveVA",
      new TGeoCompositeShape("(shVAbox + shVAHbox:trVAH1 + shVAHbox:trVAH2 + shVAC)-shVACh"),
      kMedSteelSh);
  voValve->AddNode(voValveVA, 1, 0);
  // Define Vacuum Hole of Valve
  TGeoTube   * shVACvach   = new TGeoTube("shVACvach", 0., 5.3, dz2VA/2.);
  TGeoVolume * voValveVAvh = new TGeoVolume("voValveVAvacuum", shVACvach, kMedVacuum);
  voValve->AddNode(voValveVAvh,1,0);
  // Also add valve Support (ernesto.calvo@pucp.edu.pe) 
  voValve->AddNode(voVS, 1, new TGeoTranslation(0.,  12.75, -dz2VA/2.-dzVS/2.));
  voValve->AddNode(voVS, 2, new TGeoTranslation(0., -12.75, -dz2VA/2.-dzVS/2.));

  // Define Volume VB (ernesto.calvo@pucp.edu.pe) 
  const Double_t dxVB = 23.5; 
  const Double_t dyVB =  5.0; 
  const Double_t dzVB =  9.4; 
  TGeoVolume * voValveVB = new TGeoVolume("voValveVB", 
      new TGeoBBox("shVBbox", dxVB/2., dyVB/2., dzVB/2.),
      kMedSteelSh);
  voValve->AddNode(voValveVB, 1, new TGeoTranslation(  0., dyVA/2. +dyVB/2. , 0));
  // Define Volume VC (ernesto.calvo@pucp.edu.pe) 
  const Double_t R1VC  =  4.5 /2.;
  const Double_t R2VC  =  8.1 /2.;
  const Double_t dy1VC  = 10.0;
  const Double_t dy2VC =  0.75;
  new TGeoTube("shVC1",      0.,   R1VC, dy1VC/2.);
  new TGeoTube("shVC2",      0.,   R2VC, dy2VC/2.);
  ( new TGeoTranslation("trVC21", 0., 0.,  dy1VC/2. - dy2VC/2.) )->RegisterYourself();
  ( new TGeoTranslation("trVC22", 0., 0., -dy1VC/2. + dy2VC/2.) )->RegisterYourself();
  TGeoVolume * voValveVC = new TGeoVolume("voValveVC",
      new TGeoCompositeShape("shVC1  + shVC2:trVC21 + shVC2:trVC22"),
      kMedSteelSh);
  voValve->AddNode(voValveVC, 1, new TGeoCombiTrans( 
        0., dyVA/2. + dyVB + dy1VC/2. , 0, Rx90) );
  // Define volume VD (ernesto.calvo@pucp.edu.pe) 
  const Double_t dxVD = 15.9;
  const Double_t dyVD = 23.0;
  const Double_t dzVD = 14.0;
  TGeoVolume * voValveVD = new TGeoVolume("voValveVD",
      new TGeoBBox("shVD", dxVD/2., dyVD/2., dzVD/2.),
      kMedSteelSh);
  voValve->AddNode(voValveVD, 1, 
      new TGeoTranslation( 1.25, dyVA/2. + dyVB + dy1VC + dyVD/2. , 0) );

  //
  // Define volume Front Bar (ernesto.calvo@pucp.edu.pe) 
  //
  const Double_t dxF  = 67.4; 
  const Double_t dyF  =  4.0; 
  const Double_t dzF  =  2.0; 
  const Double_t R1FC =  8.1;
  const Double_t R2FC = 11.5;
  const Double_t dxFA  = (dxF-R1FC)/2.; 
 
  new TGeoBBox("shFA", dxFA/2., dyF/2., dzF/2.);
  new TGeoTube("shFC",      0.,   R2FC, dzF/2.);
  new TGeoTube("shFCH",     0.,   R1FC, 2.*dzF);
  ( new TGeoTranslation("trFA1",  R1FC +dxFA/2., 0., 0.) )->RegisterYourself();
  ( new TGeoTranslation("trFA2", -R1FC -dxFA/2., 0., 0.) )->RegisterYourself();
  TGeoVolume * voFrontBar = new TGeoVolume("voFrontBar",
      new TGeoCompositeShape("shFA:trFA1 + shFA:trFA2 + (shFC - shFCH)"),
      kMedSteelSh);
  // Make Lateral Bars
  const Double_t kdzLatBar = 22.9;
  TGeoVolume * voLatBar = new TGeoVolume("voLatBar", 
    new TGeoTube("shLatBar",0., dyF/2., kdzLatBar/2.), kMedSteelSh);
  //
  // Define Compensator Magnet coils
  //
  dz = 12.5;
  Ro = 15.0;
  Ri = 1.9;
  new TGeoTubeSeg("shCoilRo", 0., Ro, dz/2., 90., 180.);
  new TGeoTubeSeg("shCoilRi", 0., Ri, dz   , 90., 185.);
  (new TGeoTranslation("trCoilRo", Ro, 10.4, 0.)) -> RegisterYourself();
  (new TGeoTranslation("trCoilRi", Ro,  9.6, 0.)) -> RegisterYourself();
  new TGeoBBox("shBoxCoilRo", 15.0/2., 10.4/2., dz/2.);
  new TGeoBBox("shBoxCoilRi",  1.9/2.,  9.6   , dz   );
  (new TGeoTranslation("trBoxCoilRo",       15.0/2., 10.4/2., 0.)) -> RegisterYourself();
  (new TGeoTranslation("trBoxCoilRi", 13.1 + 1.9/2.,  0.    , 0.)) -> RegisterYourself();
  new TGeoBBox("shBoxCoil0", 10./2., 30., dz);
  (new TGeoTranslation("trBoxCoil0",  14.6 + 10./2.,  0.    , 0.)) -> RegisterYourself();
  strSh  = "";
  strSh += "(shCoilRo:trCoilRo + shBoxCoilRo:trBoxCoilRo) - ";
  strSh += "(shCoilRi:trCoilRi + shBoxCoilRi:trBoxCoilRi +"  ;
  strSh += " shBoxCoil0:trBoxCoil0 )"  ;
  TGeoCompositeShape * shCoil = new TGeoCompositeShape("shCoil0", strSh);
  TGeoVolume * voCoil = new TGeoVolume("voCoil", shCoil, kMedSteelSh);

  // 
  // ALUMINIUM PLATES 
  //

  // Shape for aluminium Plate separating cavern and LHC tunnel
  const Double_t dAlWallThick = 0.5; // thickness of aluminium plates (cm)
  //
  // RB24/26 Tunnel Floor 
  R   = 220.;
  // h   = 140.;
  // phi = TMath::ACos(h / r);
  // xl  = r * TMath::Sin(phi);
  // dr  = 1600.;
  // dh  = dr * TMath::Cos(phi);
  // dl  = dr * TMath::Sin(phi);

  new TGeoTube("shWallBase",    0.,    R, dAlWallThick*0.5);  // base shape for shWallBigPlate
  new TGeoBBox("shWallCutBot",270., 110., dAlWallThick    );  // to be substracted from base
  Double_t origin[3] = {+120.+1.5+3.25, 0., 0.};
  new TGeoBBox("shWallBundleHole", 3.25,  7.5, dAlWallThick, origin);  // to be substracted from base
  // Translation for cutting circular and square hole in the plates
  (new TGeoTranslation("trAntiBeamAxis",   -70.,               40., 0.)) -> RegisterYourself();
  (new TGeoTranslation("trHUWAT3",         -70., -110. - 140. +40., 0.)) -> RegisterYourself();

  //
  //  Wall Big Aluminium Plate with Squared Hole 
  //
  const Double_t dSqrHoleSide = 33.0; // Side
  new TGeoBBox("shWallSqrHole", dSqrHoleSide*0.5, dSqrHoleSide*0.5, dAlWallThick);
  strSh  = ""; 
  strSh += "shWallBase:trAntiBeamAxis - ";
  strSh += " ( shWallCutBot:trHUWAT3" ;
  strSh += " + shWallSqrHole";
  strSh += " + shWallBundleHole)";
  TGeoVolume* voWallBigPlate = new TGeoVolume("voWallBigPlate", 
    new TGeoCompositeShape("shWallBigPlate", strSh), kMedAlu );
  //
  // Wall Squared Aluminium Plate 
  //
  const Double_t dCircHoleRad = 9.5; // Radius
  new TGeoTube("shCircHole", 0., dCircHoleRad, dAlWallThick);
  // Make holes for bars
  new TGeoTube("shRodHole", 0.,   1.7/2.,  2*dAlWallThick);
  ( new TGeoTranslation("trWallRod1",  12.5, -12.75, 0.)) -> RegisterYourself();
  ( new TGeoTranslation("trWallRod2",  12.5,  12.75, 0.)) -> RegisterYourself();
  ( new TGeoTranslation("trWallRod3", -12.5, -12.75, 0.)) -> RegisterYourself();
  ( new TGeoTranslation("trWallRod4", -12.5,  12.75, 0.)) -> RegisterYourself();
  new TGeoBBox("shWallSqrPlateBase", dSqrHoleSide*0.5 + 5.0, dSqrHoleSide*0.5 + 5.0, dAlWallThick/2.);
  strSh  = " ((((";
  strSh += " ( shWallSqrPlateBase - shCircHole )"; 
  strSh += " - shRodHole:trWallRod1)"  ;
  strSh += " - shRodHole:trWallRod2)"  ;
  strSh += " - shRodHole:trWallRod3)"  ;
  strSh += " - shRodHole:trWallRod4)"  ;
  TGeoVolume* voWallSqrPlate = new TGeoVolume("HUWAT_AlWall02", 
    new TGeoCompositeShape("shWallSqrPlate", strSh ), kMedAlu);
  // ==========================================================================
  //
  // Define Mother Vacuum volume of VMAOI  (need shIonPumpVBo)
  //
  // ==========================================================================
  const Double_t kdzMoFlange   =  2.0;
  const Double_t kdzMoBellow   = 15.6;
  const Double_t kdzTTube = 11.5; // Bellow starts here
  const Double_t kziTTube  =  1.0; // Ion Pum Tube starts here
  new TGeoTube( "shMoFlange", 0., 15.2/2., kdzMoFlange/2.0);
  new TGeoTube( "shMoBellow", 0., 11.4/2., kdzMoBellow/2.0);
  (new TGeoTranslation("trMoFlange1", 0., 0.,        0.5*kdzMoFlange     )) -> RegisterYourself();
  (new TGeoTranslation("trMoFlange2", 0., 0., 28.0 - 0.5*kdzMoFlange     )) -> RegisterYourself();
  (new TGeoTranslation("trMoBellow" , 0., 0., 28.0 - 0.5*kdzMoBellow     )) -> RegisterYourself();
  (new TGeoTranslation("trMoTTube"  , 0., 0., kziTTube +  0.5*kdzTTube   )) -> RegisterYourself();
  
  TGeoCompositeShape * shMoVMAOI =  
      new TGeoCompositeShape("shMoVMAOI", 
      "shMoFlange:trMoFlange1 + shMoBellow:trMoBellow + shMoFlange:trMoFlange2 + shIonPumpVBo:trMoTTube");
  TGeoVolume * voMoVMAOI = new TGeoVolume("voMoVMAOI", shMoVMAOI, kMedVacuum);
  //
  TGeoVolume * voVBU9mm  = new TGeoVolume("voVBU9mm" , shVBU9mm , kMedAlu);
  TGeoVolume * voVBU26mm = new TGeoVolume("voVBU26mm", shVBU26mm, kMedAlu);
  TGeoVolume * voVBUcent = new TGeoVolume("voVBUcent", shVBUcent, kMedAlu);
  voVSR      ->SetLineColor(kViolet+6);
  voIonPumpVB->SetLineColor(kViolet+6);
  voVSRflange->SetLineColor(kViolet+6);
  voVSRcontD ->SetLineColor(kViolet+6);
  voVSRcontF ->SetLineColor(kViolet+6);
  voVBUrotFlg->SetLineColor(kViolet+6); 
  voVBUrotFlg->SetLineColor(kViolet+6); 
  voVBUflg   ->SetLineColor(kViolet+6); 
  voVBU9mm   ->SetLineColor(kViolet+6);
  voVBU26mm  ->SetLineColor(kViolet+6);
  voVBUcent  ->SetLineColor(kViolet+6);
  //         
  voMoVMAOI->AddNode(voVSR      ,1, new TGeoTranslation(0., 0., 7.1 - 0.45));
  voMoVMAOI->AddNode(voIonPumpVB,1, new TGeoTranslation(0., 0., 1 + 11.5/2.));
  voMoVMAOI->AddNode(voVSRflange,1, new TGeoTranslation(0.,0.,0.));
  voMoVMAOI->AddNode(voVSRcontD ,1, new TGeoTranslation(0.,0.,28.));
  voMoVMAOI->AddNode(voVSRcontF ,1, new TGeoTranslation(0.,0.,28.));
  z = 1.0 + 11.5;
  voMoVMAOI->AddNode( voVBU9mm   , 1, new TGeoTranslation(0.,0., z) );
  voMoVMAOI->AddNode( voVBU26mm  , 1, new TGeoTranslation(0.,0., z + 11.9) );
  voMoVMAOI->AddNode( voVBUcent  , 1, new TGeoTranslation(0.,0., z) );
  voMoVMAOI->AddNode( voVBUrotFlg, 1, new TGeoCombiTrans (0.,0.,1.31, Ry180) );
  voMoVMAOI->AddNode( voVBUrotFlg, 2, new TGeoTranslation(0.,0.,28. - 1.31) );
  voMoVMAOI->AddNode( voVBUflg   , 1, new TGeoTranslation(0.,0.,0.) );
  voMoVMAOI->AddNode( voVBUflg   , 2, new TGeoCombiTrans (0.,0.,28., Ry180) );
  // ==========================================================================
  //
  // AD Support structure by Pieter Ijzerman
  // ecalvovi@cern.ch
  // ==========================================================================
  nvertices=0;
  // Cover plate_______________________________________________________________
  TGeoXtru * shADcoverplate = new TGeoXtru(2);
  shADcoverplate->SetNameTitle("shADcoverplate","shADcoverplate");
  Double_t y1[] = {  0.0, 18.50, 18.50, 22.50, 22.50, 18.50, 18.50, 22.50, 22.50, 18.50, 18.50,   .00 ,  .00, 15.25, 15.25,  .00 }; 
  Double_t x1[] = {  0.0,   .00,  5.15,  5.15, 17.15, 17.15, 24.25, 24.25, 36.25, 36.25, 41.40, 41.40 ,35.70, 35.70,  5.70, 5.70 }; 
  nvertices = sizeof(x1)/sizeof(Double_t);
  shADcoverplate->DefinePolygon(nvertices,x1,y1);
  shADcoverplate->DefineSection(0, -0.1, -20.7, 0.0, 1.0); // Z position, offset and scale for first section
  shADcoverplate->DefineSection(1,  0.1, -20.7, 0.0, 1.0); // -''- secons section

  // Horizontal side___________________________________________________________
  TGeoXtru * shADhorizontalside = new TGeoXtru(2);
  shADhorizontalside->SetNameTitle("shADhorizontalside","shADhorizontalside");
  Double_t x2[] = {  0.0,  .00, 4.80, 4.80, 7.20, 7.20, 12.00, 12.00 };
  Double_t y2[] = {  0.0, 5.66, 5.66, 1.16, 1.16, 5.66,  5.66,   .00 };
  nvertices = sizeof(x2)/sizeof(Double_t);
  shADhorizontalside->DefinePolygon(nvertices,x2,y2);
  shADhorizontalside->DefineSection(0, -0.4, -6.0, 0.0, 1.0); // Z position, offset and scale for first section
  shADhorizontalside->DefineSection(1, +0.4, -6.0, 0.0, 1.0); // -''- secons section

  TGeoBBox * shADsidebox = new TGeoBBox("shADsidebox", 0.4, 18.55/2., 5.66/2.);
  TGeoVolume * voADsidebox = new TGeoVolume("voADsidebox", shADsidebox, kMedAlu);


  TGeoVolume * voADhorizontalside = new TGeoVolume("voADhorizontalside", shADhorizontalside, kMedAlu);
  TGeoVolume * voADcoverplate = new TGeoVolume("voADcoverplate", shADcoverplate, kMedAlu);
  //
  TGeoVolume *voADsupport = new TGeoVolumeAssembly("voADsupport"); 
  voADsupport->AddNode(voADcoverplate,  1, new TGeoTranslation( 0., 0., -5.66/2.-0.23));
  voADsupport->AddNode(voADcoverplate,  2, new TGeoTranslation( 0., 0., +5.66/2.+0.23));
  voADsupport->AddNode(voADhorizontalside,  1, new TGeoCombiTrans( -6.0 - 7.1/2., 22.5-0.4, -5.66/2., Rx90));
  voADsupport->AddNode(voADhorizontalside,  2, new TGeoCombiTrans( +6.0 + 7.1/2., 22.5-0.4, -5.66/2., Rx90));
  voADsupport->AddNode(voADsidebox,  1, new TGeoTranslation( -20.7 +0.4, 18.55/2., 0.));
  voADsupport->AddNode(voADsidebox,  2, new TGeoTranslation( +20.7 -0.4, 18.55/2., 0.));
  // Add Color
  voADcoverplate     ->SetLineColor(kGray + 1);
  voADhorizontalside ->SetLineColor(kGray + 1);
  voADsidebox        ->SetLineColor(kGray + 1);
  voADsupport->SetLineColor(kGray+1);

  // ==========================================================================
  //
  // Define ADA
  //
  // ==========================================================================

  TGeoVolume *ad = new TGeoVolumeAssembly("AD");
  
  // Get medium for ADA
  TGeoMedium * medADASci        = gGeoManager->GetMedium("AD_BC404"); // AD Scin.
  // TGeoMedium * medADALG      = gGeoManager->GetMedium("AD_PMMA");  // lightGuide
  // TGeoMedium * medADAPMGlass = gGeoManager->GetMedium("AD_Glass"); // Glass for Aluminium simulation
  // TGeoMedium * medADAPMAlum  = gGeoManager->GetMedium("AD_Alum");
  
  // Get Medium for ADC 
  TGeoMedium * medADCSci     = gGeoManager->GetMedium("AD_BC404");
  // TGeoMedium * medADCLG      = gGeoManager->GetMedium("AD_PMMA");
  // TGeoMedium * medADCPMGlass = gGeoManager->GetMedium("AD_Glass");
  // TGeoMedium * medADCPMAlum  = gGeoManager->GetMedium("AD_Alum");
  
  // ADA Scintillator Pad 
  const Double_t kADACellSideY = 21.6;
  const Double_t kADACellSideX = 18.1;
  // ADC Scintillator Pad 
  const Double_t kADCCellSideY = 21.6;
  const Double_t kADCCellSideX = 18.1;
  // WLS bar          :  0.40 cm ( 4.0 mm )
  // Wrapping         :  0.20 cm ( 2.0 mm )
  // Aluminnized Mylar:  0.01 cm ( 0.1 mm )
  // Fishing line     :  0.04 cm ( 0.4 mm )
  // total shift on X :  0.65 cm
  // total shift on Y :  0.21 cm
  const Double_t kShiftX       =  0.54;
  const Double_t kShiftY       =  0.10;
  const Double_t kADACelldz    =  2.54;
  const Double_t kADCCelldz    =  2.54;
  const Double_t kADABeamPipeR =  6.20; // Radius of beam pipe hole for ADA (Diameter  12.4 cm)
  const Double_t kADCBeamPipeR =  3.70; // Radius of beam pipe hole for ADC (Diameter   7.4 cm)
  const Int_t    kColorADA     = kGreen;
  const Int_t    kColorADC     = kGreen;
  Double_t X = kShiftX + kADACellSideX * 0.5;
  Double_t Y = kShiftY + kADACellSideY * 0.5;
  Double_t WLS_dx =  0.4;
  Double_t WLS_dz =  2.5;
  Double_t WLS_SideA_Long_dy  = 24.20; // 24.2;
  Double_t WLS_SideC_Long_dy  = 24.20; // 24.2;
  Double_t WLS_SideA_Short_dy = 18.20; // 18.41; 
  Double_t WLS_SideC_Short_dy = 20.70; // 20.91; 
  // Creating ADA WLS bars_____________________________________________________
  TGeoVolume * vADA_WLS_s = new TGeoVolume( "ADAWLSshort", 
      new TGeoBBox( "shADAWLSbarShort" , WLS_dx/2.0, WLS_SideA_Short_dy/2.0, WLS_dz/2.0),
      medADASci);      
  TGeoVolume * vADA_WLS_l = new TGeoVolume( "ADAWLSlong" , 
      new TGeoBBox( "shADAWLSbarLong"  , WLS_dx/2.0, WLS_SideA_Long_dy /2.0, WLS_dz/2.0),
      medADASci);      
  vADA_WLS_l->SetLineColor( kRed );
  vADA_WLS_s->SetLineColor( kRed );
  // Creating ADC WLS bars_____________________________________________________
  TGeoVolume * vADC_WLS_s = new TGeoVolume( "ADCWLSshort", 
      new TGeoBBox( "shADCWLSbarShort" , WLS_dx/2.0, WLS_SideC_Short_dy/2.0, WLS_dz/2.0),
      medADCSci);      
  TGeoVolume * vADC_WLS_l = new TGeoVolume( "ADCWLSlong" , 
      new TGeoBBox( "shADCWLSbarLong"  , WLS_dx/2.0, WLS_SideC_Long_dy /2.0, WLS_dz/2.0),
      medADCSci);      
  vADC_WLS_l->SetLineColor(kRed);
  vADC_WLS_s->SetLineColor(kRed);
  // Make ADA scintillator pad_________________________________________________
  new TGeoBBox( "shADAbox" , kADACellSideX/2.0, kADACellSideY/2.0, kADACelldz/2.0 );
  new TGeoTube( "shADAHole",               0. , kADABeamPipeR    , kADACelldz     );
  ( new TGeoTranslation("trADAbox", X, Y, 0.)) -> RegisterYourself();
  // 
  TGeoVolume * vADA1 = new TGeoVolume( "ADApad", 
    new TGeoCompositeShape("shADApad", "shADAbox:trADAbox-shADAHole"), medADASci );      
  vADA1->SetLineColor( kColorADA );
  
  TGeoVolume *secADA  = new TGeoVolumeAssembly( "ADAsec" ); 
  // Add PAD
  Double_t fX_ADA_WLS_s = 0.1 + WLS_dx/2.0;
  Double_t fX_ADA_WLS_l = kShiftX + WLS_dx/2.0 + kADACellSideX + 0.04;
  secADA->AddNode( vADA1, 1, 0); 
  secADA->AddNode( vADA_WLS_s, 1, new TGeoTranslation(fX_ADA_WLS_s, kADABeamPipeR + WLS_SideA_Short_dy/2.0, 0.0) ); 
  secADA->AddNode( vADA_WLS_l, 1, new TGeoTranslation(fX_ADA_WLS_l,        kShiftY + WLS_SideA_Long_dy/2.0, 0.0) ); 

  /// Assembling ADA adding 4 sectors                                       //  Sectors
  TGeoVolume *vADAarray = new TGeoVolumeAssembly( "ADA" );                  //        ^ y
  vADAarray->AddNode( secADA, 1 );                                          //        |   
  vADAarray->AddNode( secADA, 2, Ry180 );                                   //   2    |   1
  vADAarray->AddNode( secADA, 3, Rz180 );                                   // --------------->  x     
  vADAarray->AddNode( secADA, 4, Rx180 );                                   //   3    |   4
  // 
  // PMT-BOX A-Side
  //
  Float_t thickbox = 0.3; // in cm
  Float_t pmbox_x = 60.;
  Float_t pmbox_y = 45.;
  Float_t pmbox_z = 10.;
  new TGeoBBox("shPMTbox1", pmbox_x/2.,           pmbox_y/2.,           pmbox_z/2.          );
  new TGeoBBox("shPMTbox2", pmbox_x/2. -thickbox, pmbox_y/2. -thickbox, pmbox_z/2. -thickbox);
  Double_t obox3[3] = {0., -2., 0.};
  new TGeoBBox("shPMTbox3", pmbox_x/2. -2*thickbox, pmbox_y/2. -thickbox, (pmbox_z-4.)/2.     , obox3);
  TGeoVolume * voPMTbox = new TGeoVolume("voPMTbox", new TGeoCompositeShape("shPMTbox","(shPMTbox1-shPMTbox2)-shPMTbox3"), kMedAlu);
  voPMTbox -> SetLineColor(kGray + 1);
  Double_t fPMTdz = 6.;
  TGeoVolume * voADApmt = new TGeoVolume("voADApmt", new TGeoTube("shADApmt", 0., 2.0, fPMTdz/2.), medADASci);
  voADApmt -> SetLineColor(kGray);
  // PMT boxes 
  printf("PMT UP y-coord: %f\n", 51.+ pmbox_y/2.); 
  printf("PMT DW y-coord: %f\n",-66 - pmbox_y/2.); 
  Double_t aX_PMT[8] = {6, -12, -12, 6, 12, -6, -6, 12};
  Double_t aY_PMT[8] = {
    51.+ pmbox_y/2.,
    51.+ pmbox_y/2.,
   -66.- pmbox_y/2.,
   -66.- pmbox_y/2.,
    51.+ pmbox_y/2.,
    51.+ pmbox_y/2.,
   -66.- pmbox_y/2.,
   -66.- pmbox_y/2.
  };
  // Layer closer to IP is shifted to the right in puicture (i.e. towards negative x axis) 
  // sector 1: Inner PMT @ x=  6cm 
  // sector 2: Inner PMT @ x=-12cm 
  // sector 3: Inner PMT @ x=-12cm
  // sector 4: Inner PMT @ x=  6cm 
  //
  // sector 1: Outer PMT @ x= 12cm 
  // sector 2: Outer PMT @ x= -6cm 
  // sector 3: Outer PMT @ x= -6cm
  // sector 4: Outer PMT @ x= 12cm 
  TGeoVolumeAssembly * voADAPMTarray = new TGeoVolumeAssembly("voADAPMTarray");
  voADAPMTarray->AddNode( voPMTbox,    1, new TGeoTranslation( 0.,  51.+ pmbox_y/2., 0.)); 
  voADAPMTarray->AddNode( voPMTbox,    2, new TGeoCombiTrans ( 0., -66 - pmbox_y/2., 0., Rz180)); 
  // PMT's
  for (Int_t i=0; i<8; i++) {
    Double_t Z1 = 0;
    if (i<4) Z1 = -kADACelldz/2. -0.23;
    else     Z1 = +kADACelldz/2. +0.23; 
    voADAPMTarray->AddNode( voADApmt, i+1, new TGeoCombiTrans(aX_PMT[i], aY_PMT[i], Z1, Rx90)); 
  }
  const Float_t kPosADA = 1696.67;  // z-center of assembly (cm) New, according to Survey by F. Klumb and E.Calvo 2015 Sept 4th.
  Double_t X1;
  Double_t Y1 = 24.3;
  Double_t Y2 = 0;
  Double_t Z1 = 0;
  for (Int_t i=0; i<8; i++) {
    Double_t sign = 1.;
    Double_t X2;
    if (i==1||i==2||i==5||i==6) sign = -1;
    X2 = aX_PMT[i];
    if (aY_PMT[i]>0) { Y1 =  24.3; Y2 = aY_PMT[i]-fPMTdz/2.; }
    else             { Y1 = -24.3; Y2 = aY_PMT[i]+fPMTdz/2.; }
    if (i<4) Z1 = -kADACelldz/2. -0.1;
    else     Z1 = +kADACelldz/2. +0.1; 
    TGeoVolume * voFiber_ADA = 0;
    fX1FiberShort[i] = sign*fX_ADA_WLS_s;
    fX1FiberLong [i] = sign*fX_ADA_WLS_l;
    fX2Fiber[i] = X2;
    fY1Fiber[i] = Y1;
    fY2Fiber[i] = Y2;
    voFiber_ADA = new TGeoVolume(Form("voFiberADAShort_%d",i), MakeFiberBundle(Form("FiberADAShort_%d",i), fX1FiberShort[i], X2, Y1, Y2, Z1), medADASci);
    voFiber_ADA -> SetLineColor(kCyan);
    voADAPMTarray->AddNode (voFiber_ADA, 1); 
    voFiber_ADA = new TGeoVolume(Form("voFiberADALong_%d",i) , MakeFiberBundle(Form("FiberADALong_%d",i) , fX1FiberLong [i], X2, Y1, Y2, Z1), medADASci);
    voFiber_ADA -> SetLineColor(kCyan);
    voADAPMTarray->AddNode (voFiber_ADA, 1); 
  }
  printf( " [ADA] %14s %14s %14s %14s %14s\n", "fX1FiberShort", "fX1FiberLong", "fX2Fiber", "fY1Fiber", "fY2Fiber" );
  for (Int_t i=0; i<8; i++) {
    printf( " [%3d] %14f %14f %14f %14f %14f\n", i, fX1FiberShort[i] , fX1FiberLong [i] , fX2Fiber[i] , fY1Fiber[i] , fY2Fiber[i] );

  }

  TGeoVolume * voADAWallPlate = new TGeoVolume("voADAWallPlate", new TGeoBBox(25., 25., 0.15), kMedSteelSh);
  voADAWallPlate -> SetLineColor(kGray+3);
  ad->AddNode ( voADAWallPlate  , 1 , new TGeoTranslation ( 0. , 0. , 1701.5                      ) ); 
  ad->AddNode ( vADAarray       , 1 , new TGeoTranslation ( 0. , 0. , kPosADA - kADACelldz/2. -0.1) ); 
  ad->AddNode ( vADAarray       , 2 , new TGeoTranslation ( 0. , 0. , kPosADA + kADACelldz/2. +0.1) ); 
  ad->AddNode ( voADsupport     , 1 , new TGeoTranslation ( 0. , 0. , kPosADA                     ) ); 
  ad->AddNode ( voADsupport     , 2 , new TGeoCombiTrans  ( 0. , 0. , kPosADA,  Rz180             ) ); 
  ad->AddNode ( voADAPMTarray   , 1 , new TGeoTranslation ( 0. , 0. , kPosADA                     ) ); 

  // ==========================================================================
  //
  // Define ADC (2014, May 4) Updated 2015, Jan 22
  //
  // ==========================================================================

  /////////////////////////////////////////////////////////////////////////////
  /// ADC in the tunnel                                                     ///
  /////////////////////////////////////////////////////////////////////////////
  new TGeoBBox( "shADCbox" , kADCCellSideX/2.0, kADCCellSideY/2.0, kADCCelldz/2.0 );
  new TGeoTube( "shADCHole",               0. , kADCBeamPipeR    , kADCCelldz     );
  X = kShiftX + kADCCellSideX * 0.5;
  Y = kShiftY + kADCCellSideY * 0.5;
  ( new TGeoTranslation("trADCbox", X, Y, 0.) ) -> RegisterYourself();
  // 
  TGeoVolume * vADCpad = new TGeoVolume( "ADCpad", 
    new TGeoCompositeShape("shADCpad", "shADCbox:trADCbox-shADCHole"), medADCSci );      
  vADCpad->SetLineColor( kColorADC );
  
  /// Creating Sector for Tunnel (Asembly:  Scintillator Pad + Light guide + PM )
  TGeoVolume *voADC  = new TGeoVolumeAssembly("ADCsec");
  // Add PAD
  voADC->AddNode( vADCpad, 1, 0);
  // Add ADC WLS Short bar
  voADC->AddNode( vADC_WLS_s, 1, 
      new TGeoTranslation( 0.1 + WLS_dx/2.0, kADCBeamPipeR + WLS_SideC_Short_dy/2.0, 0.0) ); 
  // Add ADC WLS Long  bar
  voADC->AddNode( vADC_WLS_l, 1, 
      new TGeoTranslation( 0.04 + WLS_dx/2.0 + kADCCellSideX + kShiftX, kShiftY + WLS_SideC_Long_dy/2.0, 0.0) ); 
  
  /// Assembling ADC adding the 4 sectors                 //  Sectors
  TGeoVolume *vADCarray = new TGeoVolumeAssembly("ADC");  //        ^ y
  vADCarray->AddNode( voADC, 1 );                         //        |   
  vADCarray->AddNode( voADC, 2, Ry180 );                  //   2    |   1
  vADCarray->AddNode( voADC, 3, Rz180 );                  // --------------->  x  
  vADCarray->AddNode( voADC, 4, Rx180 );                  //   3    |   4
                                                          //        |
                                                                             

  // ==========================================================================
  //
  // Add ADC to AD volume
  //
  // Note to future maintainers: 
  // In previous AliRoot versions the position z = -1900.75 corresponded 
  // to the end of the YSAA3_CC_BLOCK (concrete block shielding just before 
  // the C-Side LHC wall). Now this has been fixed to agree with reality. 
  // The YSAA3_CC_BLOCK starts at 1800.75 and ends at 1880.75 cm.
  //
  // Ernesto Calvo and Alberto Gago.
  // - ecalvovi@cern.ch
  // - agago@pucp.edu.pe
  //
  // ==========================================================================
  
  // *  ad -> AddNode(vADCarray , 1, new TGeoTranslation(0., 0., kZendADC2 + kADCCelldz/2.)); // Tunnel
  // *  ad -> AddNode(vADCarray , 2, new TGeoTranslation(0., 0., kZbegADC1 - kADCCelldz/2.)); // Tunnel
  // *  const Float_t kZbegADC1 = -kZbegFrontBar-2.
  // *  const Float_t kZendADC2 = -1959.0;            // (ecalvovi@cern.ch) 
  
  switch (fADCPosition ) {
    case kADCInTunnel:
      {
        // const Float_t kZbegADC1 = -kZbegFrontBar-2.;  // (ecalvovi@cern.ch) 
        // const Float_t kZendADC2 = -1959.0;            // (ecalvovi@cern.ch) 
        // ad -> AddNode(vADCarray , 1, new TGeoTranslation(0., 0., kZendADC2 + kADCCelldz/2.)); // Tunnel
        // ad -> AddNode(vADCarray , 2, new TGeoTranslation(0., 0., kZbegADC1 - kADCCelldz/2.)); // Tunnel
        const Float_t kPosADC = -kZbegFrontBar-2.-3.0-0.3;  // 3.0 = (5.6 + 0.2 + 0.2)/2. // (ecalvovi@cern.ch) 
        printf("CreateAD: kPosADC=%8.2f\n", kPosADC);
        // ad -> AddNode(vADCarray,   1, new TGeoTranslation(0., 0., kPosADC - kADCCelldz/2. - 0.23)); // Tunnel // ADC1
        // ad -> AddNode(vADCarray,   2, new TGeoTranslation(0., 0., kPosADC + kADCCelldz/2. + 0.23)); // Tunnel // ADC2
        // ad -> AddNode(voADsupport, 3, new TGeoTranslation(0., 0., kPosADC));
        // ad -> AddNode(voADsupport, 4, new TGeoCombiTrans (0., 0., kPosADC, Rz180));
        break;
      }
    case kADCInCavern:
      {
        printf("FATAL: vADCInCavern is now obsolete!");
        exit(1);
        break;
      }
    case kADCInBoth:
      {
        printf("FATAL: vADCInBoth   is now obsolete!");
        exit(1);
        break;
      }
  }


  // ==========================================================================
  // 
  // Add structure volumes to vADCstruct and vADAstruct volume assemblies
  //
  // ==========================================================================
  // 
  // Add some colors:
  //
  voWallBigPlate -> SetLineColor(kGray);
  voWallSqrPlate -> SetLineColor(kGray);
  voSaa3Rod      -> SetLineColor(kGray);
  voSaa3EndPlate -> SetLineColor(kGray);
  voFrontBar     -> SetLineColor(kGray);
  voLatBar       -> SetLineColor(kGray);
  voCoil         -> SetLineColor(kOrange+7);
  voValveVA      -> SetLineColor(kAzure+1);
  voValveVB      -> SetLineColor(kAzure+1);
  voValveVC      -> SetLineColor(kAzure+1);
  voValveVD      -> SetLineColor(kAzure+1);

  TGeoVolumeAssembly * vADCstruct = new TGeoVolumeAssembly("voADCStruct");
  vADCstruct->AddNode(voSaa3EndPlate, 1, new TGeoTranslation( 0., 0., kZendAbs + 1.95/2.));
  z = kZwall;
  vADCstruct->AddNode(voWallBigPlate, 1, new TGeoTranslation(0., 0., z - 0.5 * dAlWallThick ));
  vADCstruct->AddNode(voWallSqrPlate, 1, new TGeoTranslation(0., 0., z + 0.5 * dAlWallThick ));
  z = kZendAbs + 1.95 + dzRodL/2.; 
  vADCstruct->AddNode(voSaa3Rod,  1, new TGeoTranslation(  12.5, -12.75, z));
  vADCstruct->AddNode(voSaa3Rod,  2, new TGeoTranslation(  12.5,  12.75, z));
  vADCstruct->AddNode(voSaa3Rod,  3, new TGeoTranslation( -12.5, -12.75, z));
  vADCstruct->AddNode(voSaa3Rod,  4, new TGeoTranslation( -12.5,  12.75, z));
  vADCstruct->AddNode(voValve,    1, new TGeoTranslation( 0., 0., zPosValve));
  vADCstruct->AddNode(voMoVMAOI,  1, new TGeoTranslation( 0., 0., kZbegVMAOI));
  vADCstruct->AddNode(voFrontBar, 1, new TGeoTranslation( 0., 0., kZbegFrontBar + dzF/2.));
  z = kZbegCoil;
  vADCstruct->AddNode(voCoil, 1, new TGeoCombiTrans(  3.6 + dz/2., 0., z, Ry90m));
  vADCstruct->AddNode(voCoil, 2, new TGeoCombiTrans(  3.6 + dz/2., 0., z, new TGeoRotation((*Ry90m)*(*Rx180))));
  vADCstruct->AddNode(voCoil, 3, new TGeoCombiTrans( -3.6 - dz/2., 0., z, Ry90m));
  vADCstruct->AddNode(voCoil, 4, new TGeoCombiTrans( -3.6 - dz/2., 0., z, new TGeoRotation((*Ry90m)*(*Rx180))));
  z = kZbegFrontBar + dzF + kdzLatBar/2.;
  vADCstruct->AddNode(voLatBar, 1, new TGeoTranslation(  31.9, 0., z));
  vADCstruct->AddNode(voLatBar, 2, new TGeoTranslation( -31.9, 0., z));
  vADCstruct->AddNode(CreatePmtBoxC(), 1); 
  //
  // Color for voADsupp17 
  voADsupp17MainProp ->SetLineColor(kGray+1);
  voADsupp17Base     ->SetLineColor(kGray+1);
  voADsuppTopBracket ->SetLineColor(kGray+1);
  // Color for voADsupp18 
  voADsupp18MainProp -> SetLineColor(kGray+1);
  // Color for voADsuppIBeam
  voADsuppIBeam -> SetLineColor(kGreen+3);
  voADsuppIBeamV -> SetLineColor(kGreen+3);
  // voADsuppIBeam -> SetLineColorAlpha(kGreen+3, 0.);
  // voADsuppIBeam -> SetFillColorAlpha(kGreen+3, 0.5);
  // voADsuppIBeam -> SetTransparency(16);
  //
  TGeoVolume * voBLMsupport = gGeoManager->MakeBox("voBLMsupport", kMedAlu, 1.5/2.0, 30.0/2.0, 6.0/2.0);
  voBLMsupport -> SetLineColor(kGray);
  //
  TGeoVolume         * voSupportZEM           = CreateSupportZEM(); 
  // TGeoVolume         * voMagnet               = CreateMagnetYoke(); 
  TGeoVolumeAssembly * voBLM                  = CreateBLM(); 
  // TGeoVolumeAssembly * voPipeOvalyzed         = CreatePipeOvalyzed(kMedSteelSh);
  TGeoVolumeAssembly * voVacuumChamberSupport = CreateVacuumChamberSupport();
  TGeoVolumeAssembly * voShield = CreateADAShielding();
  //
  TGeoVolumeAssembly * vADAstruct = new TGeoVolumeAssembly("voADAStruct");
  vADAstruct->AddNode(voADACablingVBar      , 1 , new TGeoTranslation(28.25 ,   0.0 , 1604.00 )  );
  vADAstruct->AddNode(voADACablingVBar      , 2 , new TGeoTranslation(98.25 ,   0.0 , 1604.00 )  );
  vADAstruct->AddNode(voADAMagnetCableArray , 1 , new TGeoTranslation(26.25 ,   0.0 , 1608.00 )  );
  vADAstruct->AddNode(voADsupp17            , 1 , new TGeoTranslation(   0. , -36.5 , 1541.   )  );
  vADAstruct->AddNode(voADsupp18            , 1 , new TGeoTranslation(   0. , -36.5 , 1612.   )  ); // As Measured
  vADAstruct->AddNode(voADsupp18            , 2 , new TGeoTranslation(   0. , -36.5 , 1306.8  )  ); // As Measured 
  // vADAstruct->AddNode(voADsupp18            , 2 , new TGeoTranslation(   0. , -36.5 , 1295.8  )  ); // As LHC Drawing
  vADAstruct->AddNode(voADsuppIBeam         , 1 , new TGeoTranslation(   0. , -36.5 , 1701.   )  );
  vADAstruct->AddNode(voADsuppIBeamV        , 1 , new TGeoCombiTrans (   0. , -36.5-9.6-(81.8-9.6)/2.0,         1295.47+9.6/2.0, Rx90   )  );
  vADAstruct->AddNode(voADsuppIBeamV        , 2 , new TGeoCombiTrans (   0. , -36.5-9.6-(81.8-9.6)/2.0, 294.0 + 1295.47+9.6/2.0, Rx90   )  );
  vADAstruct->AddNode(voVacuumChamberSupport, 1 , new TGeoTranslation(   0. ,   0.0 , 1075.+125.        ));
  vADAstruct->AddNode(voVacuumChamberSupport, 2 , new TGeoCombiTrans (   0. ,   0.0 , 1075.-125. , Ry180));
  vADAstruct->AddNode(voShield              , 1 , new TGeoTranslation(   0. ,   0.0 , 1665.5  )  );
  // vADAstruct->AddNode(voPipeOvalyzed        , 1 , new TGeoTranslation(   0. ,   0.0 , 1075.47 )  );
  vADAstruct->AddNode(voBLM                 , 1 , new TGeoTranslation(  12.5,   0.0 , 1459.80 )  );
  vADAstruct->AddNode(voBLMsupport          , 1 , new TGeoTranslation(  12.5, -21.5 , 1459.80 )  );
  // vADAstruct->AddNode(voMagnet              , 1 , new TGeoTranslation(   0. ,   0.0 , 1075.00 )  );
  vADAstruct->AddNode(voSupportZEM          , 1 , new TGeoTranslation(   0. ,   0.0 ,  804.50 )  );

  if (fADCstruct) {
    ad->AddNode(vADCstruct,1, Ry180);
    ad->AddNode(vADAstruct,1       );
  }

  //
  // Add Everything to ALICE
  //
  // TGeoVolume *alice = gGeoManager->GetVolume("ALIC");
  alice->AddNode(ad, 1);
  

  // gGeoManager->DefaultColors();
  // // gGeoManager->CloseGeometry();
  // // gGeoManager->SetVisLevel(10);
  // // gGeoManager->SetVisOption(0);
  // // alice->Draw("ogl");
  /// C - Side
  return 0;
  TGLViewer * fViewer = (TGLViewer *) gPad->GetViewer3D();
  printf("TGLSAViewer: %p", fViewer);
  // from https://root.cern.ch/phpBB3/viewtopic.php?f=3&t=8143
  Double_t zoom      = 3.98;
  Double_t fov       = 2.1; // field of view
  Double_t dolly     =   30;  // distance from center
  Double_t center[3] = {0.,0.,-1670.};
  Double_t vRotate   =  65.;
  Double_t hRotate   = -10.;

  TGLViewer::ECameraType camera = TGLViewer::kCameraPerspXOZ;
  fViewer->SetCurrentCamera(camera);
  fViewer->CurrentCamera().SetExternalCenter(1);
  fViewer->CurrentCamera().SetCenterVecWarp(0., 0., 1075);
  fViewer->SetPerspectiveCamera (camera, fov, dolly, center, hRotate * TMath::DegToRad(), vRotate * TMath::DegToRad());
  /// A - Side
  ////////// // You can call gPad->GetViewer3D() and cast this to TGLViewer (or TGLSAViewer)
  ////////// // TGLSAViewer * fViewer = (TGLSAViewer *) gPad->GetViewer3D();
  ////////// TGLViewer * fViewer = (TGLViewer *) gPad->GetViewer3D();
  ////////// printf("TGLSAViewer: %p", fViewer);
  ////////// fViewer->CurrentCamera().SetExternalCenter(1);
  ////////// fViewer->CurrentCamera().SetCenterVecWarp(0., 0., 1690);
  ////////// // from https://root.cern.ch/phpBB3/viewtopic.php?f=3&t=8143
  ////////// Double_t zoom      = 0.98;
  ////////// Double_t fov       = 3.5; // field of view
  ////////// Double_t dolly     = 3;  // distance from center
  ////////// Double_t center[3] = {0,0,1690};
  ////////// Double_t vRotate   = 75.;
  ////////// Double_t hRotate   = 0;

  ////////// TGLViewer::ECameraType camera = TGLViewer::kCameraPerspXOZ;
  ////////// fViewer->SetCurrentCamera(camera);
  ////////// fViewer->SetPerspectiveCamera (camera, fov, dolly, center, hRotate * TMath::DegToRad(), vRotate * TMath::DegToRad());
  printf("<=== AliADv1::CreateAD(): ver=[Feb 3st, 2015]; contact=[ecalvovi@cern.ch]\n");
  return 0; 
  //
  // TGLViewer::ECameraType camera = TGLViewer::kCameraOrthoXOZ;
  // TGLViewer *v = (TGLViewer *) gPad->GetViewer3D();
  // v->SetCurrentCamera(camera);
  // v->CurrentCamera().SetExternalCenter(kTRUE);

  // Double_t zoom      = 0.98;
  // Double_t dolly     = 4;
  // Double_t center[3] = {0,0,1690};
  // Double_t hRotate   = 0.01;
  // Double_t vRotate   = 0.01;

  // v->SetOrthoCamera(camera, zoom, dolly, center, hRotate, vRotate);
  // v->DoDraw();
}
*/

//_________________________________________________________


TGeoVolumeAssembly * CreateVacuumChamberSupport()
{
  // Drawings used in this section:
  // LHC Drawing: LHCVC2U_0009
  // LHC Drawing: LHCVC2U_0018
  // LHC Drawing: LHCVC2U_0019
  // LHC Drawing: LHCVC2U_0014
  // LHC Drawing: LHCVC2U_0016
  // Originally  made by: Sergio Best.
  // Edited and Reviewed: Ernesto Calvo.
  // Support should start +- 125 cm from magnet center
  //
  TGeoVolumeAssembly * voVacuumChamberSupport = new TGeoVolumeAssembly("voVacuumChamberSupport");
  Float_t zpos = 19.;
  /////////////////////////////////////////////////////////////////////////////
  //
  // LHC Drawing: LHCVC2U_0013
  //
  TGeoXtru * sh_LHCVC2U_0013_a = NULL;
  TGeoBBox * sh_LHCVC2U_0013_b = NULL;
  // TGeoVolume * vol_LHCVC2U_0013 = new TGeoVolume("vol_LHCVC2U_0013", sh_LHCVC2U_0013, kMedSteelSh);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // sh_LHCVC2U_0013_a:
  sh_LHCVC2U_0013_a = new TGeoXtru(2);
  Double_t xpoints[] = { 0. , 0.95 , 22.45 , 24.45 , 24.45 , 22.45 , 0.95 , 0.};
  Double_t ypoints[] = { 0. ,   0. , 8.25  , 8.25  , 9.75  , 9.75  , 3.   , 3.};
  Int_t nvertices = sizeof(xpoints)/sizeof(Double_t);
  sh_LHCVC2U_0013_a -> DefinePolygon(nvertices,xpoints,ypoints);
  sh_LHCVC2U_0013_a -> DefineSection(0, -1.5, 0., 0., 1.0); // index, Z position, offset (x,y) and scale for first section
  sh_LHCVC2U_0013_a -> DefineSection(1,  0.0, 0., 0., 1.0); // idem, second section
  sh_LHCVC2U_0013_a -> SetName("sh_LHCVC2U_0013_a");

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // sh_LHCVC2U_0013_b:
  sh_LHCVC2U_0013_b = new TGeoBBox("sh_LHCVC2U_0013_b", 2.75/2., 3.0/2., 0.95/2.);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // a+b: sh_LHCVC2U_0013
  ( new TGeoCombiTrans ("tr_a" , 0.               , 0.     , 0.      ,  Ry90m )) -> RegisterYourself();
  ( new TGeoTranslation("tr_b" , +1.5 + 2.75/2.   , 3.0/2. , 0.95/2.          )) -> RegisterYourself();
  TGeoVolume * vol_LHCVC2U_0013   = new TGeoVolume("vol_LHCVC2U_0013", 
      new TGeoCompositeShape("sh_LHCVC2U_0013", "sh_LHCVC2U_0013_a:tr_a + sh_LHCVC2U_0013_b:tr_b"),
      kMedSteelSh);
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // a+b: sh_LHCVC2U_0010
  ( new TGeoCombiTrans ("tr_a10" , -1.5           , 0.     , 0.      ,  Ry90m )) -> RegisterYourself();
  ( new TGeoTranslation("tr_b10" , -1.5 - 2.75/2. , 3.0/2. , 0.95/2.          )) -> RegisterYourself();
  TGeoVolume * vol_LHCVC2U_0010   = new TGeoVolume("vol_LHCVC2U_0010", 
      new TGeoCompositeShape("sh_LHCVC2U_0010", "sh_LHCVC2U_0013_a:tr_a10 + sh_LHCVC2U_0013_b:tr_b10"),
      kMedSteelSh);
  
  /////////////////////////////////////////////////////////////////////////////
  //
  // LHC Drawing: LHCVC2U_0009
  //
  new TGeoBBox("shADclampSup_0" , 5.1/2.  , 7./2. , 1.5/2.);
  new TGeoBBox("shADclampSup_1" , 2./2.   , 7./2. , 9./2. );
  new TGeoPara("shADclampSup_2" , 2./2.   , 7./2. , 20./2.   , 0. , +180./TMath::Pi()*atan(12./200.) , 0.);
  new TGeoBBox("shADclampSup_3" , 2./2.   , 7./2. , 10./2.);
  new TGeoBBox("shADclampSup_4" , 5.35/2. , 7./2. , 1.5/2.);

  TGeoRotation * RZ180      = new TGeoRotation("RZ180",   180., 0.,   0.);

  (new TGeoTranslation("trADclampSup_0" , 3.10/2. +3.2 , 0. ,    39. - 1.5/2. ))-> RegisterYourself();
  (new TGeoTranslation("trADclampSup_1" ,         +2.2 , 0. ,    30. + 9.0/2. ))-> RegisterYourself();
  (new TGeoTranslation("trADclampSup_2" ,         +1.6 , 0. ,    10. +20.0/2.  ))-> RegisterYourself();
  (new TGeoTranslation("trADclampSup_3" ,         +1.0 , 0. ,        +10.0/2. ))-> RegisterYourself();
  (new TGeoTranslation("trADclampSup_4" , 5.35/2. +0.0 , 0. ,        + 1.5/2. ))-> RegisterYourself();

  TGeoVolume* vol_LHCVC2U_0009 = new TGeoVolume("vol_LHCVC2U_0009", 
      // new TGeoCompositeShape("cshape0", "shADclampSup_1:trADclampSup_1 + shADclampSup_2:trADclampSup_2 + shADclampSup_3:trADclampSup_3 + shADclampSup_4:trADclampSup_4"),
      new TGeoCompositeShape("cshape0", "shADclampSup_0:trADclampSup_0 + shADclampSup_1:trADclampSup_1 + shADclampSup_2:trADclampSup_2 + shADclampSup_3:trADclampSup_3 + shADclampSup_4:trADclampSup_4"),
      kMedSteelSh);
  

  /////////////////////////////////////////////////////////////////////////////
  //
  // LHC Drawing: LHCVC2U_0018
  //

  new TGeoBBox("shADspacer_0",  4.5/2., 5.5/2.,  8.92/2.);
  new TGeoBBox("shADspacer_1",  4.5/2., 1.5/2.,  2.08/2.);
  new TGeoBBox("shADspacer_2",  4.1/2., 7.0/2.,  2.5 /2.);
  new TGeoBBox("shADspacer_3",  0.2/2., 5.5/2.,  2.5 /2.);

  (new TGeoTranslation("trADspacerL_0" , 0.                 ,  0.            , 0.             )) ->RegisterYourself();
  (new TGeoTranslation("trADspacerL_1" , 0.                 ,  (1.5-5.5)/2.  , 11./2.         )) ->RegisterYourself();
  (new TGeoTranslation("trADspacerL_2" , -4.5/2.-0.2-4.1/2. ,  (7.0-5.5)/2.  , -(8.92-2.5)/2. )) ->RegisterYourself();
  (new TGeoTranslation("trADspacerL_3" , -4.5/2.-0.1        ,            0.  , -(8.92-2.5)/2. )) ->RegisterYourself();


  TGeoVolume* vol_LHCVC2U_0018 = new TGeoVolume("vol_LHCVC2U_0018", 
      // new TGeoCompositeShape("cshape2", "shADspacer_0:trADspacerL_0+shADspacer_1:trADspacerL_1+shADspacer_2:trADspacerL_2"), kMedSteelSh);
      new TGeoCompositeShape("cshape2", "shADspacer_0:trADspacerL_0+shADspacer_1:trADspacerL_1+shADspacer_2:trADspacerL_2+shADspacer_3:trADspacerL_3"), kMedSteelSh);

  /////////////////////////////////////////////////////////////////////////////
  //
  // LHC Drawing: LHCVC2U_0019
  //
  (new TGeoTranslation("trADspacerR_0" , 0.                 ,  0.            , 0.             )) ->RegisterYourself();
  (new TGeoTranslation("trADspacerR_1" , 0.                 ,  (1.5-5.5)/2.  , 11./2.         )) ->RegisterYourself();
  (new TGeoTranslation("trADspacerR_2" , +4.5/2.+0.2+4.1/2. ,  (7.0-5.5)/2.  , -(8.92-2.5)/2. )) ->RegisterYourself();
  (new TGeoTranslation("trADspacerR_3" , +4.5/2.+0.1        ,            0.  , -(8.92-2.5)/2. )) ->RegisterYourself();

  // (new TGeoTranslation("trADspacerR_1",  0.,(1.5+5.5)/2.5,11./2. ))->RegisterYourself();
  // (new TGeoTranslation("trADspacerR_2",  8.8/2.,-(7.-5.5)/2.,-(8.92-2.5)/2. ))->RegisterYourself();

  TGeoVolume* vol_LHCVC2U_0019 = new TGeoVolume("vol_LHCVC2U_0019", 
      new TGeoCompositeShape("cshape3", "shADspacer_0:trADspacerR_0+shADspacer_1:trADspacerR_1+shADspacer_2:trADspacerR_2+shADspacer_3:trADspacerR_3"), kMedSteelSh);
      // new TGeoCompositeShape("cshape3", "shADspacer_0:trADspacerR_0+shADspacer_1:trADspacerR_1+shADspacer_2:trADspacerR_2"), kMedSteelSh);

  /////////////////////////////////////////////////////////////////////////////
  //
  // LHC Drawing: LHCVC2U_0014
  //
  /////////////////////////////////////////////////////////////////////////////

  Float_t h2_inf=(7.5*sqrt(3.)-7.6)/2.;
  Float_t h1=7.6-h2_inf;
  Float_t trp_m=20.-(2.*h1)/sqrt(3.);

  TGeoRotation * RW90      = new TGeoRotation("RW90",   0., -90.,   0.);

  new TGeoBBox("shADInf_0" , 30./2.   , 2./2.    , 4./2.  ); 
  new TGeoTrd1("shADInf_1" , 20./2.   , trp_m/2. , 2./2.  ,  h1/2.     ); 
  new TGeoTrd1("shADInf_2" , trp_m/2. , 5./2.    , 2./2.  ,  h2_inf/2. ); 
  new TGeoTube("shADInf_3" , 0.       , 7.62     , 10./2. ); 

  new TGeoBBox("shADInf_4" , 25./2.   , 3./2.    , 10./2. ); 

  (new TGeoTranslation("trADInf_0" , 0. , 0. , 0.                          )) ->RegisterYourself();
  (new TGeoTranslation("trADInf_1" , 0. , 0. , 2. + (h1/2.)                )) ->RegisterYourself();
  (new TGeoTranslation("trADInf_2" , 0. , 0. , 2. + h1+(h2_inf/2.)         )) ->RegisterYourself();
  (new TGeoCombiTrans ("trADInf_3" , 0. , 0. , 0.0                  , RW90 )) ->RegisterYourself();
  (new TGeoTranslation("trADInf_4" , 0. , 0. , -(10./2.)+0.1               )) ->RegisterYourself();

  TGeoVolume* vol_LHCVC2U_0014 = new TGeoVolume("vol_LHCVC2U_0014", 
      new TGeoCompositeShape("shADInf_0:trADInf_0+shADInf_1:trADInf_1+shADInf_2:trADInf_2-shADInf_3:trADInf_3-shADInf_4:trADInf_4"), kMedSteelSh);


  /////////////////////////////////////////////////////////////////////////////
  //
  // LHC Drawing: LHCVC2U_0016
  //
  /////////////////////////////////////////////////////////////////////////////

  TGeoRotation * RW90m      = new TGeoRotation("RW90",   0., 90.,   0.);

  new TGeoBBox("shADSup_0" , 24.8/2.  , 2./2.    , 1.9/2.              ); 
  new TGeoTrd1("shADSup_1" , 20./2.   , trp_m/2. , 2./2.  ,  h1/2.     ); 
  new TGeoTrd1("shADSup_2" , trp_m/2. , 5./2.    , 2./2.  ,  h2_inf/2. ); 
  new TGeoTube("shADSup_3" , 0.       , 7.62     , 10./2.              ); 

  (new TGeoTranslation("trADSup_0" , 0. , 0. , 0.1 + 1.9/2.           )) ->RegisterYourself();
  (new TGeoTranslation("trADSup_1" , 0. , 0. , (4.+h1)/2.             )) ->RegisterYourself();
  (new TGeoTranslation("trADSup_2" , 0. , 0. , (4.+2.*h1+h2_inf)/2.   )) ->RegisterYourself();
  (new TGeoCombiTrans ("trADSup_3" , 0. , 0. , 0. , RW90      )) ->RegisterYourself();

  TGeoVolume * vol_LHCVC2U_0016 = new TGeoVolume("vol_LHCVC2U_0016", 
      new TGeoCompositeShape("(shADSup_0:trADSup_0+shADSup_1:trADSup_1+shADSup_2:trADSup_2)-shADSup_3:trADSup_3"), kMedSteelSh);

  vol_LHCVC2U_0009->SetLineColor(kCyan     );  // Clamp Support
  vol_LHCVC2U_0010->SetLineColor(kOrange-3 );  // Left  Stiffener
  vol_LHCVC2U_0013->SetLineColor(kOrange-3 );  // Right Stiffener
  vol_LHCVC2U_0016->SetLineColor(kViolet+6 );  // Top      semi-circular support
  vol_LHCVC2U_0014->SetLineColor(kViolet+8 );  // Inferior semi-circular support
  vol_LHCVC2U_0019->SetLineColor(kGreen    );  // Right Spacer
  vol_LHCVC2U_0018->SetLineColor(kSpring-5 );  // Left Spacer

  voVacuumChamberSupport->AddNode ( vol_LHCVC2U_0009 , 1 , new TGeoCombiTrans  (  -3. , 0.00        , 0.              , RZ180 )); 
  voVacuumChamberSupport->AddNode ( vol_LHCVC2U_0009 , 2 , new TGeoCombiTrans  (  +3. , 0.00        , 0.              , NULL  )); 
  voVacuumChamberSupport->AddNode ( vol_LHCVC2U_0010 , 1 , new TGeoCombiTrans  ( -4.1 , -14. + 0.75 , 0.              , NULL  )); 
  voVacuumChamberSupport->AddNode ( vol_LHCVC2U_0013 , 1 , new TGeoCombiTrans  ( +4.1 , -14. + 0.75 , 0.              , NULL  )); 
  voVacuumChamberSupport->AddNode ( vol_LHCVC2U_0016 , 1 , new TGeoCombiTrans  ( 0.   , 0.          , zpos+10.+9.+11. , RW90  )); 
  voVacuumChamberSupport->AddNode ( vol_LHCVC2U_0014 , 1 , new TGeoCombiTrans  ( 0.   , 0.          , zpos+10.+9.+11. , RW90m )); 
  voVacuumChamberSupport->AddNode ( vol_LHCVC2U_0019 , 1 , new TGeoTranslation ( -15.+4.5/2. , -0.75        , zpos+20.+8.92/2.    )); 
  voVacuumChamberSupport->AddNode ( vol_LHCVC2U_0018 , 1 , new TGeoTranslation (  15.-4.5/2. , -0.75        , zpos+20.+8.92/2.    )); 

  return voVacuumChamberSupport;

}

TGeoVolumeAssembly * CreateADAShielding()
{
  TGeoVolumeAssembly * voShieldADA = new TGeoVolumeAssembly("voShieldADA");
  TGeoVolumeAssembly * voShieldTop = new TGeoVolumeAssembly("voShieldADA_Top");
  TGeoVolumeAssembly * voShieldBot = new TGeoVolumeAssembly("voShieldADA_Bot");
  // This section made by Ernesto Calvo (PUCP)
  // Base Plate (x1):
  // dx: 45.5 cm
  // dz: 33.0 cm
  // dy:  0.5 cm
  //
  // Support bars along x axis (x3):
  // 50x50x48 cm thickness: ??
  //

  // Shield Cage:
  // Height  : 25.0 cm
  // Inner dx: 40.0 cm
  // Inner dz: 30.0 cm
  //
  // F1,2: 50.0 x 25.0 x  0.5 cm
  // E1,2:  0.5 x 25.0 x 30.0 cm
  // H1,2:  4.0 x 25.0 x  0.5 cm
  TGeoVolumeAssembly * voShieldCage = new TGeoVolumeAssembly("voADAShieldCage");
  TGeoVolume * volF = gGeoManager -> MakeBox("volAD_Sh_F", kMedSteelSh, 50.0/2., 25.0/2.,  0.5/2.); 
  TGeoVolume * volE = gGeoManager -> MakeBox("volAD_Sh_E", kMedSteelSh,  0.5/2., 25.0/2., 30.0/2.); 
  TGeoVolume * volH = gGeoManager -> MakeBox("volAD_Sh_H", kMedSteelSh,  4.0/2., 25.0/2.,  0.5/2.); 
  TGeoVolume * volS = gGeoManager -> MakeBox("volAD_Sh_S", kMedSteelSh, 40.0/2., 30.0/2., 30.0/2.); 
  // C1  : 47.0 x  5.0 x  5.0 cm
  new TGeoBBox("shRecTubC_outer" , 47./2. , 5./2.     , 5./2.     ); 
  new TGeoBBox("shRecTubC_inner" , 47.    , 5./2.-0.3 , 5./2.-0.3 ); 
  new TGeoBBox("shRecTubD_outer" , 46./2. , 5./2.     , 5./2.     ); 
  new TGeoBBox("shRecTubD_inner" , 46.    , 5./2.-0.3 , 5./2.-0.3 ); 
  // B1  : 4.0 x  6.0 x  13.0 cm
  new TGeoBBox("shRecTubB1_outer" , 4./2.     , 6./2.     , 13./2. ); 
  new TGeoBBox("shRecTubB1_inner" , 4./2.-0.3 , 6./2.-0.3 , 13.    ); 
  // B2  : 4.0 x  6.0 x  15.0 cm
  new TGeoBBox("shRecTubB2_outer" , 4./2.     , 6./2.     , 15./2. ); 
  new TGeoBBox("shRecTubB2_inner" , 4./2.-0.3 , 6./2.-0.3 , 15.    ); 
  // A   : 5.0 x  5.0 x  57.0 cm
  new TGeoBBox("shRecTubA_outer" , 5./2.      , 57./2.    , 5./2.      ); 
  new TGeoBBox("shRecTubA_inner" , 5./2.-0.3  , 57.       , 5./2. -0.3 ); 
  //
  TGeoVolume * volC  = new TGeoVolume("volAD_Sh_C"  , new TGeoCompositeShape("shRecTubC"  , "shRecTubC_outer-shRecTubC_inner"  ) , kMedSteelSh ); 
  TGeoVolume * volD  = new TGeoVolume("volAD_Sh_D"  , new TGeoCompositeShape("shRecTubD"  , "shRecTubD_outer-shRecTubD_inner"  ) , kMedSteelSh ); 
  TGeoVolume * volB1 = new TGeoVolume("volAD_Sh_B1" , new TGeoCompositeShape("shRecTubB1" , "shRecTubB1_outer-shRecTubB1_inner") , kMedSteelSh ); 
  TGeoVolume * volB2 = new TGeoVolume("volAD_Sh_B2" , new TGeoCompositeShape("shRecTubB2" , "shRecTubB2_outer-shRecTubB2_inner") , kMedSteelSh ); 
  TGeoVolume * volA  = new TGeoVolume("volAD_Sh_A"  , new TGeoCompositeShape("shRecTubA"  , "shRecTubA_outer-shRecTubA_inner"  ) , kMedSteelSh ); 
  // IBEAMS:
  const Double_t IBeamLength = 200.;
  TGeoVolume * volIBeamTop2 = MakeVolIBeam("volAD_Sh_IBeamTop2", kMedSteelSh, 18., 18., 0.8, 1.4, IBeamLength);
  TGeoVolume * volIBeamTop  = MakeVolIBeam("volAD_Sh_IBeamTop" , kMedSteelSh, 10., 10., 0.6, 1.0, 58.0);
  TGeoVolume * volIBeamBot  = MakeVolIBeam("volAD_Sh_IBeamBot" , kMedSteelSh, 10., 10., 0.6, 1.0, 33.0);
  volIBeamBot -> SetLineColor(kGreen+3);
  volIBeamTop -> SetLineColor(kGreen+3);
  volIBeamTop2-> SetLineColor(kGreen-3);

  TGeoVolume * volPlateTop = gGeoManager -> MakeBox("volAD_Sh_PlateTop", kMedSteelSh, 45.5/2.,  0.5/2., 33.0/2.); 
  TGeoVolume * volPlateBot = gGeoManager -> MakeBox("volAD_Sh_PlateBot", kMedSteelSh, 50.0/2.,  0.5/2., 33.0/2.); 
  volF -> SetLineColor(kOrange);
  volE -> SetLineColor(kOrange-8);
  volH -> SetLineColor(kOrange-8);
  volS -> SetLineColor(kGray+3);
  volC -> SetLineColor(kGreen-3);
  volD -> SetLineColor(kGreen-3);
  volB1-> SetLineColor(kRed-10);
  volB2-> SetLineColor(kRed-10);
  volA -> SetLineColor(kTeal);
  volPlateTop -> SetLineColor(kViolet+7);
  volPlateBot -> SetLineColor(kViolet+7);

  voShieldCage->AddNode(volF , 1 , new TGeoTranslation(         0.0 , 25.0/2. , -15. -0.25 )); 
  voShieldCage->AddNode(volF , 2 , new TGeoTranslation(         0.0 , 25.0/2. , +15. +0.25 )); 
  voShieldCage->AddNode(volE , 1 , new TGeoTranslation( -20.0 -0.25 , 25.0/2. , 0.00       )); 
  voShieldCage->AddNode(volE , 2 , new TGeoTranslation( +20.0 +0.25 , 25.0/2. , 0.00       )); 
  voShieldCage->AddNode(volH , 1 , new TGeoTranslation( +22.0 +0.50 , 25.0/2. , -15.0+0.25 )); 
  voShieldCage->AddNode(volH , 2 , new TGeoTranslation( -22.0 -0.50 , 25.0/2. , -15.0+0.25 )); 
  voShieldCage->AddNode(volH , 3 , new TGeoTranslation( +22.0 +0.50 , 25.0/2. , +15.0-0.25 )); 
  voShieldCage->AddNode(volH , 4 , new TGeoTranslation( -22.0 -0.50 , 25.0/2. , +15.0-0.25 )); 
  voShieldCage->AddNode(volS , 1 , new TGeoTranslation(         0.0 , 30.0/2. , 0.00       )); 

  voShieldTop -> AddNode(voShieldCage , 1 , new TGeoTranslation(  0.   , 0.    , 0.         )); 
  voShieldTop -> AddNode(volPlateTop  , 1 , new TGeoTranslation( -2.25 , -0.25 , +1.0       )); 
  voShieldTop -> AddNode(volC         , 1 , new TGeoTranslation( -2.50 , -3.00 , -12.5      )); 
  voShieldTop -> AddNode(volC         , 2 , new TGeoTranslation( -2.50 , -3.00 , -12.5+25.0 )); 
  voShieldTop -> AddNode(volD         , 1 , new TGeoTranslation( -2.50 , -3.00 ,  0.0       )); 
  voShieldTop -> AddNode(volB1        , 1 , new TGeoTranslation(+23.00 , -3.00 , -9.0       )); 
  voShieldTop -> AddNode(volB1        , 2 , new TGeoTranslation(-28.00 , -3.00 , -9.0       )); 
  voShieldTop -> AddNode(volB2        , 1 , new TGeoTranslation(+23.00 , -3.00 ,+10.0       )); 
  voShieldTop -> AddNode(volB2        , 2 , new TGeoTranslation(-28.00 , -3.00 ,+10.0       )); 
  voShieldTop -> AddNode(volA         , 1 , new TGeoTranslation(+23.00 ,+20.50 ,  0.0       )); 
  voShieldTop -> AddNode(volA         , 2 , new TGeoTranslation(-28.00 ,+20.50 ,  0.0       )); 
  voShieldTop -> AddNode(volA         , 3 , new TGeoTranslation(+23.00 ,+20.50 ,+20.0       )); 
  voShieldTop -> AddNode(volA         , 4 , new TGeoTranslation(-28.00 ,+20.50 ,+20.0       )); 
  voShieldTop -> AddNode(volIBeamTop  , 1 , new TGeoCombiTrans ( -2.50 ,+60.00 ,  0.0, Ry90 )); 
  voShieldTop -> AddNode(volIBeamTop  , 2 , new TGeoCombiTrans ( -2.50 ,+60.00 ,+20.0, Ry90 )); 
  Double_t trz = +(2.25 + IBeamLength/2. -15.5 );
  voShieldTop -> AddNode(volIBeamTop2 , 1 , new TGeoTranslation( +0.00 ,+50.00 ,trz         )); 

  //
  voShieldBot -> AddNode(voShieldCage , 1 , new TGeoCombiTrans (0.    , 0.    , 0.  , Ry180 )); 
  voShieldBot -> AddNode(volPlateBot  , 1 , new TGeoTranslation(0.    , -0.25 , -1.0        )); 
  voShieldBot -> AddNode(volIBeamBot  , 1 , new TGeoTranslation(-20.0 , -0.50 , -1.0        )); 
  voShieldBot -> AddNode(volIBeamBot  , 2 , new TGeoTranslation(0.    , -0.50 , -1.0        )); 
  voShieldBot -> AddNode(volIBeamBot  , 3 , new TGeoTranslation(+20.0 , -0.50 , -1.0        )); 
  //
  voShieldADA -> AddNode(voShieldTop  , 1 , new TGeoCombiTrans (0.    ,   82.0 +2.0 -30.0, 0.  , Ry180 ));
  voShieldADA -> AddNode(voShieldBot  , 1 , new TGeoTranslation(0.    , -124.6 -2.0      , 0.          ));
  return voShieldADA;
}

TGeoVolume * MakeVolIBeam(const char * volname, const TGeoMedium * mat, const Double_t x, const Double_t y, const Double_t dx, const Double_t dy, const Double_t dz)
{
  TGeoXtru * shIBeam = new TGeoXtru(2);
  TString shname = "sh";
  shname += volname;
  shIBeam->SetNameTitle(shname.Data(), shname.Data());
  Double_t H, H2, W, W2;
  W = x/2.; W2 = dx/2.; 
  H = y/2.; H2 = dy/1.; 
  Double_t xIBeam[] = { -W , -W   , -W2  , -W2  , -W   , -W , W , W    , W2   , W2   , W    , W  }; 
  Double_t yIBeam[] = { -H , H2-H , H2-H , H-H2 , H-H2 , H  , H , H-H2 , H-H2 , H2-H , H2-H , -H }; 
  Int_t nvertices = sizeof(xIBeam)/sizeof(Double_t);
  shIBeam -> DefinePolygon(nvertices,xIBeam,yIBeam);
  shIBeam -> DefineSection(0, -dz/2., 0., -H, 1.0); // index, Z position, offset (x,y) and scale for first section
  shIBeam -> DefineSection(1,  dz/2., 0., -H, 1.0); // idem, second section
  
  return new TGeoVolume(volname, shIBeam, mat);
}

TGeoVolume * MakeEllipticalSlice(Double_t x, Double_t y, Double_t dz, Double_t thk, Int_t iplane, const TGeoMedium * mat)
{
  TString siplane   = TString::Format("%03d", iplane);
  TString volname   = "voTubeSliceElliptical_"; volname   += siplane;
  TString shapename = "shTubeSliceElliptical_"; shapename += siplane;
  TString shapestr  = TString::Format("shTubeSliceEltuO_%s-shTubeSliceEltuI_%s", siplane.Data(), siplane.Data()); 

  new TGeoEltu(Form("shTubeSliceEltuO_%s" , siplane.Data()) , x     , y     , dz/2.);
  new TGeoEltu(Form("shTubeSliceEltuI_%s" , siplane.Data()) , x-thk , y-thk , dz/1.);


  TGeoVolume * vol =  new TGeoVolume(volname.Data(), 
                                     new TGeoCompositeShape(shapename.Data(), shapestr.Data()), 
                                     mat);
  vol->SetLineColor(kViolet+1+(iplane%2));
  return vol;
}


TGeoVolumeAssembly * CreatePipeOvalyzed(const TGeoMedium * mat)
{
  TGeoVolumeAssembly * voPipeOvalyzed = new TGeoVolumeAssembly("voRB242Pipe");

  // Create Central Part

  // TODO: ***OVERLAP***
  // const Double_t dx        =   6.75/2.0; // <-- According to https://edms.cern.ch/ui/file/909866/0/lhcvc2u_0036-v0.pdf
  const Double_t ax        =   6.00/2.0; // According to LHCVC2U_0008 ( there is a 0.75 cm discrepancy)
  const Double_t by        =   9.50/2.0;
  const Double_t dL        = 270.00/2.0;
  const Double_t thickness =   0.40/2.0;

  new TGeoEltu("shPipeOvalyzedOuter" , ax           , by           , dL          ); 
  new TGeoEltu("shPipeOvalyzedInner" , ax-thickness , by-thickness , dL+thickness); 
  TGeoCompositeShape * shPipeOvalyzedCentral = new TGeoCompositeShape("shPipeOvalyzedCentral", "shPipeOvalyzedOuter-shPipeOvalyzedInner");
  TGeoVolume * voPipeOvalyzedCentral = new TGeoVolume("voPipeOvalyzedCentral", shPipeOvalyzedCentral, mat);

  // Create Transition Elliptic->Circular: 


  // Formula: y(z) = a + (r-a)*(z/L)
  TGeoVolumeAssembly * voPipeTransition = new TGeoVolumeAssembly("voPipeTransition");

  const Double_t       Fdz =                   2.0/1.0 ;  // z-lenght of the Beam-Pipe Flange
  const Double_t        FR =                  15.2/2.0 ;  // Outher Radius of Beam-Pipe Flange
  const Double_t         R =                   8.4/2.0 ;  // Outher Radius of Beam-Pipe
  const Double_t         L =                  20.0/1.0 ;  // length of elliptical -> circular transition in cm
  const Int_t      nslices =                 100       ; 
  const Double_t        dz = L / Double_t(nslices)     ; 
  const Double_t         a = ax;
  const Double_t         b = by;

  for (Int_t iplane=0; iplane<nslices; iplane++) 
  {
    Double_t zslice = dz * Double_t(iplane);

    Double_t y = b + (R-b)*(zslice/L);
    Double_t x = a + (R-a)*(zslice/L);
    voPipeTransition->AddNode(MakeEllipticalSlice(x, y, dz, thickness, iplane, mat), 1, new TGeoTranslation(0., 0., zslice+dz/2.0));
  }

  // PIPE END
  new TGeoTube("shEndTube"  , R-thickness,  R,   L/2.0);
  new TGeoTube("shEndFlange", R-thickness, FR, Fdz/2.0);
  (new TGeoTranslation("trEndF", 0., 0., L/2.0-Fdz/2.0))->RegisterYourself();

  TGeoVolume * voRB242PipeEnd = new TGeoVolume("voRB242PipeEnd", new TGeoCompositeShape("shRB242PipeEnd", "shEndTube+shEndFlange:trEndF"), mat);
  voRB242PipeEnd -> SetLineColor(kSpring-8);

  
  /*
  TGeoXtru * shPipeTransition = new TGeoXtru(2);

  const Int_t nvertices =   360.;
  const Double_t dtheta = ( 360./ Double_t(nvertices));
  Double_t X[nvertices] = {   0. };
  Double_t Y[nvertices] = {   0. };

  for (Int_t iplane=0; iplane<nplanes; iplane++) 
  {
    for (Int_t it=0; it<nvertices; it++)
    {
      Double_t angle = Double_t(it) * dtheta;
      X[it] = a * TMath::Cos(angle); 
      Y[it] = b * TMath::Sin(angle); 
    }
    shPipeTransition -> DefinePolygon(nvertices,X,Y);
    shPipeTransition -> DefineSection(0, -dz/2., 0., -H, 1.0); // index, Z position, offset (x,y) and scale for first section
  }
  */
  voPipeOvalyzedCentral -> SetLineColor(kSpring+5);
  voPipeOvalyzed->AddNode(voPipeOvalyzedCentral, 1, new TGeoTranslation(0., 0., 0.           ));
  voPipeOvalyzed->AddNode(voPipeTransition     , 1, new TGeoTranslation(0., 0.,   dL         ));
  voPipeOvalyzed->AddNode(voPipeTransition     , 1, new TGeoCombiTrans (0., 0., -(dL ), Ry180));
  voPipeOvalyzed->AddNode(voRB242PipeEnd       , 1, new TGeoTranslation(0., 0., + dL + (1.5*L)        ));
  voPipeOvalyzed->AddNode(voRB242PipeEnd       , 1, new TGeoCombiTrans (0., 0., - dL - (1.5*L), Ry180 ));
  return voPipeOvalyzed;
}

TGeoVolumeAssembly * CreateBLM()
{
  TGeoVolumeAssembly * voBLM = new TGeoVolumeAssembly("voBLM");
  TGeoMedium * kMedAlu       = gGeoManager->GetMedium("AD_Alum");   // Aluminium 
  TGeoMedium * kMedVacuum    = gGeoManager->GetMedium("AD_VA_C0");  // Stainless Steel 
  TGeoMedium * kMedSteelSh   = gGeoManager->GetMedium("AD_ST_C0");  

  const Double_t cm  =              1.0; // 1 mm = 0.1 cm
  const Double_t mm  =              0.1; // 1 mm = 0.1 cm
  const Double_t deg = TMath::Pi()/180.; // 
// Cilinder
  Double_t SourTubeInnerRadius =   0.   * mm             ; // full tube
  Double_t SourTubeOuterRadius =  45.   * mm             ; //
  //   Double_t SourTubeLenght = 615.   * mm             ; //
  Double_t SourTubeLenght      = 495.   * mm             ; // without electronicsbox
  Double_t SourTubeHalfLenght  =   0.5  * SourTubeLenght ; // half size
  Double_t SourTubexpos        =   0.                    ; // -61. * mm * sin(GetRotAngle()) ; //0. ; 
  Double_t SourTubezpos        =   0.                    ; //  61. * mm * cos(GetRotAngle()) ; //0. ; 

  TGeoVolume * solidSourTube = gGeoManager->MakeTube( "voSourChamber", kMedVacuum, SourTubeInnerRadius, SourTubeOuterRadius, SourTubeHalfLenght);
  //
  // outer tube	inox
  //
  Double_t OuterTubeInnerRadius =  42.45*mm;
  Double_t OuterTubeOuterRadius =  44.45*mm;
  Double_t OuterTubeLenght      = 585.00*mm; // As Measured!
  // Double_t OuterTubeLenght      = 615.00*mm; // As Measured!
  // Double_t OuterTubeLenght      = 483.  *mm; // Original
  Double_t OuterTubeHalfLenght  =   0.5 *OuterTubeLenght;	// half size
  //   Double_t OuterTubezpos = -60.*mm;
  Double_t OuterTubezpos = 0.*mm;		// without electronicsbox

  TGeoVolume * voOuterTube = gGeoManager->MakeTube("voOuterTube", kMedSteelSh, OuterTubeInnerRadius, OuterTubeOuterRadius, OuterTubeHalfLenght);
// Create gas
// Create Electrode With holes
// Place the 61 electrodes
//
  // Electronic box (NO)
  // new TGeoBBox("shElectronicBoxOuter", 0.5*4.5*cm, 0.5*6.0*cm, 0.5*6.0*cm); //its size
  // new TGeoBBox("shElectronicBoxInner", 0.5*3.7*cm, 0.5*5.2*cm, 0.5*5.2*cm); //its size
  // TGeoCompositeShape * shElectronicBox = new TGeoCompositeShape("shElectronicBox", "shElectronicBoxOuter-shElectronicBoxInner");
  // TGeoVolume * voElectronicBox = new TGeoVolume("voElectronicBox", shElectronicBox, kMedAlu);
  //
  // top plate diff inox (part 6)
  //
  const Double_t PlateInnerRadius = 0.  * mm           ; 
  const Double_t PlateWidth       = 5.  * mm           ; 
  const Double_t PlateHalfLenght  = 0.5 * PlateWidth ; // half size
  const Double_t InsulHalfLenght  = 0.5 * 15. * mm ; // half size

  // FIXME: define Aluminim Oxide Al2O3
  TGeoMedium * kMedAl2O3 = kMedSteelSh;

  TGeoVolume * voEndPlateBLM = gGeoManager->MakeTube("voEndPlateBLM",	kMedSteelSh, 0.0, OuterTubeOuterRadius, PlateHalfLenght);
  TGeoVolume * voInsPlateBLM = gGeoManager->MakeTube("voInsPlateBLM",	kMedAl2O3  , 0.0, OuterTubeInnerRadius, InsulHalfLenght);



  //
  // electrodes
  //
  const Double_t ElInnerRadius = 0.   * mm;
  const Double_t ElOuterRad    = 37.5 * mm;
  const Double_t Elwidt        = 0.5  * mm;
  const Double_t Elhalfwidth   = 0.5  * Elwidt;
  const Double_t TubfixRadius  = 31.5 * mm;
  

  TGeoTube * solidelectrode = new TGeoTube("shElect",	ElInnerRadius, ElOuterRad, Elhalfwidth); 
      
  //
  // transformation
  //
  Double_t xhole = 0.;
  Double_t yhole = 0.;

  xhole = TubfixRadius * TMath::Cos( 30.*deg );
  yhole = TubfixRadius * TMath::Sin( 30.*deg );

  TGeoCompositeShape * shElectrode;
  //
  // hole in electrodes
  //
  Double_t holeInnerRadius = 0.0 * mm;
  Double_t holeOuterRad    = 3.0 * mm;
  Double_t holelenght      = 1.0 * mm;
  Double_t holehalflenght  = 0.5 * holelenght;
  TGeoTube * solidHole = new TGeoTube("hole", holeInnerRadius, holeOuterRad, holehalflenght);
  //
  for(Int_t rep=0; rep<6; rep++)
  {
    xhole = TubfixRadius * TMath::Cos( 30.*deg+rep*60.*deg );
    yhole = TubfixRadius * TMath::Sin( 30.*deg+rep*60.*deg );
    (new TGeoTranslation(Form("TrElectrodeHole%d", rep), xhole,yhole,0.)) -> RegisterYourself();
  }

  shElectrode = new TGeoCompositeShape(
      "shElectrode", 
      "shElect-(hole+"
            "hole:TrElectrodeHole0+"
            "hole:TrElectrodeHole1+"
            "hole:TrElectrodeHole2+"
            "hole:TrElectrodeHole3+"
            "hole:TrElectrodeHole4+"
            "hole:TrElectrodeHole5)");


  TGeoVolume * voElectrode = new TGeoVolume("electrode", shElectrode, kMedAlu);
  //
  // Rods
  //
  const Double_t SDGasLenght = 356.5*mm;
  const Double_t RodLength = SDGasLenght;
  const Double_t RodHalfLength = 0.5*RodLength;
  const Double_t deltaPhi = 30.* deg;

  TGeoVolume * voInnerRodBLM = gGeoManager->MakeTube("voInnerRodBLM",	kMedSteelSh, 0.0  , 0.2, RodHalfLength);
  TGeoVolume * voOuterRodBLM = gGeoManager->MakeTube("voOuterRodBLM",	kMedAlu    , 0.205, 0.3, RodHalfLength);
  voInnerRodBLM -> SetLineColor(kGray+1);
  voOuterRodBLM -> SetLineColor(kGray  );

  TGeoVolumeAssembly * voRodBLM = new TGeoVolumeAssembly("voRodBLM");
  voRodBLM->AddNode(voInnerRodBLM, 1, new TGeoTranslation(0., 0., 0.));
  voRodBLM->AddNode(voOuterRodBLM, 1, new TGeoTranslation(0., 0., 0.));


  for (Int_t Rodrep=0; Rodrep<6; Rodrep++)
  {
    Double_t xRod = TubfixRadius * TMath::Cos( deltaPhi+(Rodrep%6)*60.*deg );
    Double_t yRod = TubfixRadius * TMath::Sin( deltaPhi+(Rodrep%6)*60.*deg );

    voBLM -> AddNode(voRodBLM, Rodrep, new TGeoTranslation(xRod, yRod, 0.));
  }


  
  // voElectronicBox -> SetLineColor(kGray     ); 
  voEndPlateBLM   -> SetLineColor(kGray + 1 ); 
  voInsPlateBLM   -> SetLineColor(kYellow-10); 
  voElectrode     -> SetLineColor(kGray     ); 
  voOuterTube     -> SetLineColor(kOrange-2 ); 
  //
  // placeing of the electrodes
  //
  Double_t zEl = 0.;
  Double_t ElParamLength = 345.*mm;

  for(Int_t Elparam=0;Elparam<61;Elparam++)
  {
    zEl = -0.5*ElParamLength + 5.75*mm*Elparam;
    voBLM -> AddNode(voElectrode, Elparam+1, new TGeoTranslation(0.,0.,zEl));
  }

  // voBLM -> AddNode(voElectronicBox , 1 , new TGeoTranslation(0. , 0. , -OuterTubeHalfLenght+6.0+0.1                      )); 
  voBLM -> AddNode(voOuterTube     , 1 , new TGeoTranslation(0. , 0. , OuterTubezpos                                     )); 
  voBLM -> AddNode(voEndPlateBLM   , 1 , new TGeoTranslation(0. , 0. , OuterTubezpos-PlateHalfLenght-OuterTubeHalfLenght )); 
  voBLM -> AddNode(voEndPlateBLM   , 2 , new TGeoTranslation(0. , 0. , OuterTubezpos+PlateHalfLenght+OuterTubeHalfLenght )); 
  voBLM -> AddNode(voInsPlateBLM   , 1 , new TGeoTranslation(0. , 0. , -InsulHalfLenght-RodHalfLength                    )); 
  voBLM -> AddNode(voInsPlateBLM   , 2 , new TGeoTranslation(0. , 0. , +InsulHalfLenght+RodHalfLength                    )); 
  return voBLM;  
}


TGeoVolume * CreateWarmModule()
{
//************************************************************************************************************************************************************

//******************************************************************************************************************************************

  TGeoVolumeAssembly * ST0521057_VAMPA = new TGeoVolumeAssembly("ST0521057_VAMPA");
  /////////////////////////////////////////////////////////////////////////////
  //
  // c:\Users\jsestak\Desktop\ST0521057_VAMPA.stp /marvin.ascencio@pucp.edu.pe
  // Original System CATIA Version 5-6 Relase 2013 SP5
  // Organizarion CERN
  // https://a360.autodesk.com/viewer/#id/dXJuOmFkc2sub2JqZWN0czpvcy5vYmplY3Q6YTM2MHZpZXdlci90MTQ5NTAxNzM5MzczNF8wNjU0NTAxMzE2NDkzMzU5M18xNDk1MDE3MzkzNzM1LnN0cA
  
  const Double_t tcks_        =   0.3;//thickness
  const Double_t tcks_tube    =   0.3;//thickness
  
  TGeoBBox * VP_0002_Box1   = new TGeoBBox("VP_0002_Box1", 13.6/2,  23.5/2, 16.0/2 );
  TGeoBBox * VP_0002_Box1_aux1   = new TGeoBBox("VP_0002_Box1_aux1", (13.6/2-2*tcks_-5.6/2)*0.5,  23.5/2+0.1, 11.2/2+0.3/2 );//X2
  TGeoBBox * VP_0002_Box1_aux2   = new TGeoBBox("VP_0002_Box1_aux2", (13.6/2-2*tcks_) ,  23.5/2+0.1,14.2/2-tcks_-11.2/2 );
  TGeoBBox * VP_0002_Box1_aux3   = new TGeoBBox("VP_0002_Box1_aux3", 5.6/2 - 2*tcks_ ,  23.5/2-2*tcks_,11.2/2-2*tcks_ );
  
  //forwar boxes
  TGeoBBox * VP_0002_Box1_f1   = new TGeoBBox("VP_0002_Box1_f1", 12.8/2, 3.9/2, 0.4/2 );
  TGeoBBox * VP_0002_Box1_f2   = new TGeoBBox("VP_0002_Box1_f2", 11.2/2, 3.1/2, 1.4/2 );
  TGeoBBox * VP_0002_Box1_f3   = new TGeoBBox("VP_0002_Box1_f3", 5.6/2, 2.1/2,  11.2/2 );
  
  
  //backward boxes
  TGeoBBox * VP_0002_Box1_b1   = new TGeoBBox("VP_0002_Box1_b1",  12.8/2, 2.6/2, 0.4/2 );
  TGeoBBox * VP_0002_Box1_b2   = new TGeoBBox("VP_0002_Box1_b2", 11.2/2,  1.8/2, 1.4/2 );
  TGeoBBox * VP_0002_Box1_b3   = new TGeoBBox("VP_0002_Box1_b3", 5.6/2, 1.5/2, 11.2/2);

  //TUBES__VP_002_1
  TGeoTube * VP_0002_Tube_001     = new TGeoTube("VP_0002_Tube_001"     , 0.    ,  1.9867  , 20.8/2 );//major_tube
  TGeoTube * VP_0002_Tube_002     = new TGeoTube("VP_0002_Tube_002"     , 0.    ,  1.7141  , 1.2/2 );//intermediate_tube
  TGeoTube * VP_0002_Tube_003     = new TGeoTube("VP_0002_Tube_003"     , 0.    ,  0.8974  , 2.1/2 );//minor tube

  //TUBES__VP_002_2
  TGeoTube * VP_0002_Tube_021     = new TGeoTube("VP_0002_Tube_021"     , (5.05257-tcks_tube)    ,  5.05257  , 11.8/2 ); //huge tube for VP_002
  TGeoTube * VP_0002_Tube_021_b     = new TGeoTube("VP_0002_Tube_021_b"     , 0.0    ,  (5.05257-tcks_tube)  , tcks_tube/2 ); 
  TGeoTube * VP_0002_Tube_021_     = new TGeoTube("VP_0002_Tube_021_"     , 0.0    ,  5.05257  , 11.8/2 ); //huge tube for VP_002 FOR CUT
  TGeoTube * VP_0002_Tube_022     = new TGeoTube("VP_0002_Tube_022"     , 5.6419-3.1564+ tcks_tube   ,  5.6419  , 1.8/2 );  //short tube with large area at VP_002
  TGeoTube * VP_0002_Tube_023     = new TGeoTube("VP_0002_Tube_023"     ,  (3.1564- tcks_tube)    ,  3.1564  , 2.3/2);  //This tube will used for cut : 1.2/2 
  TGeoTube * VP_0002_Tube_023_c     = new TGeoTube("VP_0002_Tube_023_c"     ,  0.0    ,  3.1564  , 2.3/2+4);  //This tube will used for cut : 1.2/2 


  //TUBE FOR CUT
/*  TGeoTube * VP_0002_Tube_001     = new TGeoTube("VP_0002_Tube_001"     , 5.05257- tcks_tube   ,  5.05257  , 20.8/2 );//major_tube*/
  
  
  (new TGeoCombiTrans("trTube", -5.10257,  0., 0.0, Ry90m)) -> RegisterYourself();
  (new TGeoCombiTrans("trTube_tb1", -5.10257+3,  0., -16.0/2 -7.3, Ry90m)) -> RegisterYourself();
  (new TGeoTranslation("trTube_tb2", 0.0, 0.0, -16.0/2 - 11.8/2)) -> RegisterYourself();
  (new TGeoTranslation("trTube_tb3", 0.0, 0.0, 0.0)) -> RegisterYourself();
//(new TGeoTranslation("trTube_cut",-5.05257+2.3/2-1.2+5.10257+tcks_tube, 0.0, -16.0/2 -7.3)) -> RegisterYourself();
  /*(new TGeoCombiTrans("trTUBE", TGeoCombiTrans(-5.0-5.6/2-0.5*(13.6-4*tcks_-5.6)*0.5, 0.0, -15.9/2+11.2/2+1.8)5257+2.3/2-1.2, 0.0, 11.8/2 -7.3,Ry90m)) -> RegisterYourself();*/
  (new TGeoTranslation("trTube_3",0.0, 0.0, 9.8+8.0+4.8/2)) -> RegisterYourself();
  (new TGeoTranslation("trTube_4",5.6/2+0.5*(13.6-4*tcks_-5.6)*0.5, 0.0, -15.9/2+11.2/2+1.8+0.3)) -> RegisterYourself();
  (new TGeoTranslation("trTube_5",-5.6/2-0.5*(13.6-4*tcks_-5.6)*0.5, 0.0, -15.9/2+11.2/2+1.8+0.3)) -> RegisterYourself();
  (new TGeoTranslation("trTube_6",0.0, 0.0, -15.9/2+11.2/2+1.8+11.2/2+0.5*(3.0-2*tcks_))) -> RegisterYourself();
  
  TGeoCompositeShape * VP_0002_Tube_023_ = new TGeoCompositeShape("VP_0002_Tube_023_", "VP_0002_Tube_023:trTube - VP_0002_Tube_021_");
  TGeoCompositeShape * Hole_Box01 = new TGeoCompositeShape("Hole_Box01", "VP_0002_Box1 - VP_0002_Box1_aux2:trTube_6 - VP_0002_Box1_aux1:trTube_4 - VP_0002_Box1_aux1:trTube_5 - VP_0002_Box1_aux3 ");
  TGeoCompositeShape * VP_0002_Tube_021_c = new TGeoCompositeShape ("VP_0002_Tube_021_c", "VP_0002_Tube_021:trTube_tb2 - VP_0002_Tube_023_c:trTube_tb1");
  
  
  TGeoVolume * VP_0002_01 = new TGeoVolume("VP_0002_01", VP_0002_Box1, kMedSteelSh);
  VP_0002_01 -> SetLineColor(kGray);
  
  TGeoVolume * VP_0002_02 = new TGeoVolume("VP_0002_02", VP_0002_Box1_f1, kMedSteelSh);
  VP_0002_02 -> SetLineColor(kGray);

  TGeoVolume * VP_0002_03 = new TGeoVolume("VP_0002_03", VP_0002_Box1_f2, kMedSteelSh);
  VP_0002_03 -> SetLineColor(kGray);

  TGeoVolume * VP_0002_04 = new TGeoVolume("VP_0002_04", VP_0002_Box1_f3, kMedSteelSh);
  VP_0002_04 -> SetLineColor(kGray);
  
  TGeoVolume * VP_0002_05 = new TGeoVolume("VP_0002_05", VP_0002_Box1_b1, kMedSteelSh);
  VP_0002_05 -> SetLineColor(kGray);
  
  TGeoVolume * VP_0002_06 = new TGeoVolume("VP_0002_06", VP_0002_Box1_b2, kMedSteelSh);
  VP_0002_06 -> SetLineColor(kGray);   
  
  TGeoVolume * VP_0002_07 = new TGeoVolume("VP_0002_07", VP_0002_Box1_b3, kMedSteelSh);
  VP_0002_07 -> SetLineColor(kGray);

  TGeoVolume * VP_0002_08 = new TGeoVolume("VP_0002_08", VP_0002_Tube_003, kMedSteelSh);
  VP_0002_08 -> SetLineColor(kGray);
  
  TGeoVolume * VP_0002_09 = new TGeoVolume("VP_0002_09", VP_0002_Tube_002, kMedSteelSh);
  VP_0002_09 -> SetLineColor(kGray);

  TGeoVolume * VP_0002_10 = new TGeoVolume("VP_0002_10", VP_0002_Tube_001, kMedSteelSh);
  VP_0002_10 -> SetLineColor(kGray);

  TGeoVolume * VP_0002_11 = new TGeoVolume("VP_0002_11", VP_0002_Tube_021_c, kMedSteelSh);
  VP_0002_11 -> SetLineColor(kGray);
  
  TGeoVolume * VP_0002_11_a = new TGeoVolume("VP_0002_11_a", VP_0002_Tube_021_b, kMedSteelSh);
  VP_0002_11_a -> SetLineColor(kGray);

  TGeoVolume * VP_0002_12 = new TGeoVolume("VP_0002_12", VP_0002_Tube_022, kMedSteelSh);
  VP_0002_12 -> SetLineColor(kGray);

  TGeoVolume * VP_0002_13 = new TGeoVolume("VP_0002_13", VP_0002_Tube_023_, kMedSteelSh);
  VP_0002_13 -> SetLineColor(kGray);

  TGeoVolume * VP_0002_14 = new TGeoVolume("VP_0002_14", VP_0002_Box1_aux1, kMedSteelSh);
  VP_0002_14 -> SetLineColor(kGray);

  TGeoVolume * VP_0002_15 = new TGeoVolume("VP_0002_15", VP_0002_Box1_aux2, kMedSteelSh);
  VP_0002_15 -> SetLineColor(kGray);

  TGeoVolume * VP_0002_16 = new TGeoVolume("VP_0002_16", Hole_Box01, kMedSteelSh);
  VP_0002_16 -> SetLineColor(kGray);  

  
  
////////////////////////////////////////////////////////////////////////////////  
//PARTBODY
  TGeoBBox * PartBody_bx01   = new TGeoBBox("PartBody_bx01", 1.5/2,  15.0/2, 16.0/2 );  
//SHORT TUBES
  TGeoTube * PartBody_t001     = new TGeoTube("PartBody_t001"     , 0.    ,   7.5694 , 2.0/2 );//need two
  TGeoTube * PartBody_t002     = new TGeoTube("PartBody_t002"     , 5.67004-3.5593+tcks_tube    ,   5.67004 , 1.8/2 );//**-->BOX
  TGeoTube * PartBody_t003     = new TGeoTube("PartBody_t003"     , 0.    ,   5.67004 , 1.7/2 );//*#
  TGeoTube * PartBody_t004     = new TGeoTube("PartBody_t004"     , 0.    ,   3.4779 , 1.3/2 );
      
//LONG TUBES
  TGeoTube * PartBody_t005     = new TGeoTube("PartBody_t005"     , 3.1816-tcks_tube    ,   3.1816 , 7.0/2 );//Cut to 5.9 with box 
  TGeoTube * PartBody_t006     = new TGeoTube("PartBody_t006"     , 3.5593-tcks_tube    ,   3.5593 , 6.7/2 );//Cut to 5.6 with circle * #
  TGeoTube * PartBody_t007     = new TGeoTube("PartBody_t007"     , 3.5593 - tcks_tube    ,   3.5593 , 6.7/2 );//Cut to 5.6 with circle **
  TGeoTube * PartBody_t008     = new TGeoTube("PartBody_t008"     , 0.    ,   1.8881 , 5.1/2 );//cut to 4.8
  
//HUGE TUBES 
  TGeoTube * PartBody_t009     = new TGeoTube("PartBody_t009"     , 0.    ,   5.1183 , 11.1/2 );
  TGeoTube * PartBody_t010     = new TGeoTube("PartBody_t010"     , 0.    ,   5.1183 , 2.3/2 );
  TGeoTube * PartBody_t011     = new TGeoTube("PartBody_t011"     , 0.    ,   5.825 , 10.4/2 );
  
//CONE
TGeoCone * PartBody_t012     = new TGeoCone("PartBody_t012", 0.4767/2, 0.0 , 5.1183 , 0.0, 5.825);
TGeoCone * PartBody_t013     = new TGeoCone("PartBody_t013", 0.4767/2, 0.0 , 5.825 , 0.0, 5.1183);
  
  (new TGeoCombiTrans("trTube_1",-13.2026, 0.0, -15.3,Ry90m)) -> RegisterYourself();//tube_small
  (new TGeoTranslation("trTube_2",-20.5709, 0.0, -16.4093)) -> RegisterYourself();//tube_small
  (new TGeoCombiTrans("trTube_31",-28.0892, 0.0, -15.3,Ry90m)) -> RegisterYourself();
  
  /*(new TGeoCombiTrans("trTUBE", TGeoCombiTrans(-5.05257+2.3/2-1.2, 0.0, 11.8/2 -7.3,Ry90m)) -> RegisterYourself();*/

  TGeoCompositeShape * PartBody_t014 = new TGeoCompositeShape("PartBody_t014", "PartBody_t007:trTube_1 - PartBody_t009:trTube_2");//remember for cut
  
  TGeoCompositeShape * PartBody_t005_c = new TGeoCompositeShape("PartBody_t005_c", "PartBody_t005:trTube_31 - PartBody_t009:trTube_2");//remember for cut
 
  TGeoVolume * PartBody_001 = new TGeoVolume("PartBody_001", PartBody_t002, kMedSteelSh);
  PartBody_001 -> SetLineColor(kGray);

  TGeoVolume * PartBody_002 = new TGeoVolume("PartBody_002", PartBody_t014, kMedSteelSh);
  PartBody_002 -> SetLineColor(kGray);  
  
  TGeoVolume * PartBody_003 = new TGeoVolume("PartBody_003", PartBody_t009, kMedSteelSh);
  PartBody_003 -> SetLineColor(kGray);

  TGeoVolume * PartBody_004 = new TGeoVolume("PartBody_004", PartBody_t001, kMedSteelSh);
  PartBody_004 -> SetLineColor(kGray);

  TGeoVolume * PartBody_005 = new TGeoVolume("PartBody_005", PartBody_t012, kMedSteelSh);
  PartBody_005 -> SetLineColor(kGray);

  TGeoVolume * PartBody_006 = new TGeoVolume("PartBody_006", PartBody_t011, kMedSteelSh);
  PartBody_006 -> SetLineColor(kGray);

  TGeoVolume * PartBody_007 = new TGeoVolume("PartBody_007", PartBody_t013, kMedSteelSh);
  PartBody_007 -> SetLineColor(kGray);

  TGeoVolume * PartBody_008 = new TGeoVolume("PartBody_008", PartBody_t010, kMedSteelSh);
  PartBody_008 -> SetLineColor(kGray);

  TGeoVolume * PartBody_009 = new TGeoVolume("PartBody_009", PartBody_t001, kMedSteelSh);
  PartBody_009 -> SetLineColor(kGray);

  TGeoVolume * PartBody_010 = new TGeoVolume("PartBody_010", PartBody_t005_c, kMedSteelSh);
  PartBody_010 -> SetLineColor(kGray);
  
  TGeoVolume * PartBody_011 = new TGeoVolume("PartBody_011", PartBody_bx01, kMedSteelSh);
  PartBody_011 -> SetLineColor(kGray);
  //*********************** 
  TGeoVolumeAssembly * RB241_WarmModule = new TGeoVolumeAssembly("RB241_WarmModule");
  RB241_WarmModule -> AddNode(VP_0002_16, 1 , new TGeoTranslation(0., 0., 0.));
  RB241_WarmModule -> AddNode(VP_0002_02, 1 , new TGeoTranslation(0., 23.5/2+3.9/2, -15.9/2+0.4));
  RB241_WarmModule -> AddNode(VP_0002_03, 1 , new TGeoTranslation(0., 23.5/2+3.1/2, -15.9/2+1.0));
  RB241_WarmModule -> AddNode(VP_0002_04, 1 , new TGeoTranslation(0., 23.5/2+2.1/2, -15.9/2+11.2/2+1.8));
  RB241_WarmModule -> AddNode(VP_0002_05, 1 , new TGeoTranslation(0., -23.5/2-2.6/2, -15.9/2+0.4));
  RB241_WarmModule -> AddNode(VP_0002_06, 1 , new TGeoTranslation(0., -23.5/2-1.8/2, -15.9/2+1.0));
  RB241_WarmModule -> AddNode(VP_0002_07, 1 , new TGeoTranslation(0., -23.5/2-1.5/2, -15.9/2+11.2/2+1.8));

  RB241_WarmModule -> AddNode(VP_0002_08, 1 , new TGeoTranslation(-11.2/2+0.8974, 23.5/2+(3.1-0.8974), -15.9/2+2.1/2+1.8));
  RB241_WarmModule -> AddNode(VP_0002_09, 1 , new TGeoTranslation(-11.2/2+0.8974, 23.5/2+(3.1-0.8974), -15.9/2+1.2/2+2.1+1.8));
  RB241_WarmModule -> AddNode(VP_0002_10, 1 , new TGeoTranslation(-11.2/2+0.8974, 23.5/2+(3.1-0.8974), -15.9/2+1.2+20.8/2+2.1+1.8));
  RB241_WarmModule -> AddNode(VP_0002_11, 1 , new TGeoTranslation(0.0, 0.0, 0.0));
  RB241_WarmModule -> AddNode(VP_0002_11_a, 1 , new TGeoTranslation(0.0, 0.0, -16.0/2 - 11.8+tcks_tube/2));
  RB241_WarmModule -> AddNode(VP_0002_12, 1 , new TGeoCombiTrans(-5.05257-1.2-1.8/2, 0.0, -16.0/2 -7.3,Ry90m));
 RB241_WarmModule -> AddNode(VP_0002_13, 1 , new TGeoTranslation(-5.05257+2.3/2-1.2+5.10257, 0.0, -16.0/2 -7.3));
  
 //******************** 
 /* RB241_WarmModule -> AddNode(VP_0002_14, 1 , new TGeoTranslation(5.6/2+0.5*(13.6-4*tcks_-5.6)*0.5, 0.0, -15.9/2+11.2/2+1.8));
  RB241_WarmModule -> AddNode(VP_0002_14, 2 , new TGeoTranslation(-5.6/2-0.5*(13.6-4*tcks_-5.6)*0.5, 0.0, -15.9/2+11.2/2+1.8));
  RB241_WarmModule -> AddNode(VP_0002_15, 1 , new TGeoTranslation(0.0, 0.0, -15.9/2+11.2/2+1.8+11.2/2+0.5*(3.0-tcks_)));*/
 /* RB241_WarmModule -> AddNode(VP_0002_16, 1 , new TGeoTranslation(0.0, 0.0, 0.0));*/
  
//PartBody
  RB241_WarmModule -> AddNode(PartBody_001, 1 , new TGeoCombiTrans(-5.05257-1.2-1.8-1.8/2, 0.0, -16.0/2 -7.3,Ry90m));
  RB241_WarmModule -> AddNode(PartBody_002, 1 , new TGeoTranslation(0.0,0.0,0.0));
  RB241_WarmModule -> AddNode(PartBody_003, 1 , new TGeoTranslation(-5.05257-1.2-1.8-1.8-6.7+1.1-5.1183, 0.0, -16.0/2 -7.3-1.1093));
  RB241_WarmModule -> AddNode(PartBody_004, 1 , new TGeoTranslation(-5.05257-1.2-1.8-1.8-6.7+1.1-5.1183, 0.0, -16.0/2 -7.3-1.1093-11.1/2-2.0/2));
  RB241_WarmModule -> AddNode(PartBody_005, 1 , new TGeoTranslation(-5.05257-1.2-1.8-1.8-6.7+1.1-5.1183, 0.0, -16.0/2 -7.3-1.1093+11.1/2+0.4767/2));  
  RB241_WarmModule -> AddNode(PartBody_006, 1 , new TGeoTranslation(-5.05257-1.2-1.8-1.8-6.7+1.1-5.1183, 0.0, -16.0/2 -7.3-1.1093+11.1/2+0.4767+10.4/2));
  RB241_WarmModule -> AddNode(PartBody_007, 1 , new TGeoTranslation(-5.05257-1.2-1.8-1.8-6.7+1.1-5.1183, 0.0, -16.0/2 -7.3-1.1093+11.1/2+0.4767+10.4+0.4767/2));
  RB241_WarmModule -> AddNode(PartBody_008, 1 , new TGeoTranslation(-5.05257-1.2-1.8-1.8-6.7+1.1-5.1183, 0.0, -16.0/2 -7.3-1.1093+11.1/2+0.4767+10.4+0.4767+2.3/2));
  RB241_WarmModule -> AddNode(PartBody_009, 1 , new TGeoTranslation(-5.05257-1.2-1.8-1.8-6.7+1.1-5.1183, 0.0, -16.0/2 -7.3-1.1093+11.1/2+0.4767+10.4+0.4767+2.3+2.0/2));
  RB241_WarmModule -> AddNode(PartBody_010, 1 , new TGeoTranslation(0.0 , 0.0 , 0.0 ));
 //RB241_WarmModule -> AddNode(PartBody_010, 1 , new TGeoCombiTrans(-5.05257 - 1.2 - 1.8 - 1.8 - 5.6 - 5.1183*2 - 2.4, 0.0, -15.3, Ry90m));
 RB241_WarmModule -> AddNode(PartBody_011, 1 , new TGeoTranslation(-5.05257 - 1.2 - 1.8 - 1.8 - 5.6 - 5.1183*2 -5.9-1.5/2, 0.0, -15.3+3.7816));
  
 //RB241_WarmModule -> AddNode(VP_0002_13, 1 , new TGeoCombiTrans(-5.05257+2.3/2-1.2, 0.0, -16.0/2 -7.3,Ry90m));       
 /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 ///
 ///
   /////////////////////////////////////////////////////////////////////////////
  //
  // Make Vacuum Suport Vab foot base plate15 in RB24/1 FP (B-Side) /marvin.ascencio@pucp.edu.pe
  // Flat Bar
  // https://edms.cern.ch/ui/file/373462/AA/lhcvhl__0007-vAA.pdf
  // https://edms.cern.ch/file/439918/AA/lhcvhl__0027-vAA.pdf
  // Foot tube: LHCVHL__0007
  // LHC Drawing: LHCVHL__0005
  Double_t paf    = TMath::Pi()/4;
  Double_t sqrtt  = TMath::Sqrt(3);
  Double_t ds     = 2* (TMath::Sin(paf));
  Double_t fxt    = -25 + (3/sqrtt - 0.9);
  Double_t fxt2   = -25 - (3/sqrtt + 3.9);
  Double_t fxtd   =  25 - (3/sqrtt - 0.9);
  Double_t fxtd2  =  25 + (3/sqrtt + 3.9);
  Double_t fyt    =  - (sqrtt*(3/sqrtt -0.9));
  Double_t fyt2   =  (sqrtt*(3/sqrtt - 0.9)) ;
  Double_t trFx[] = {fxt , -25 , fxt2, -25. };
  Double_t trFy[] = {0.  , fyt, 0, fyt2 };
  Double_t trFx2[] = {fxtd , 25. , fxtd2, 25. };
  Double_t trFy2[] = {0.  , fyt, 0, fyt2 };
  TGeoBBox * ShADFlatBar   = new TGeoBBox("ShADFlatBar"   , 50/2  ,  35/2 , 2/2   ); 
  TGeoBBox * ShADFlatBarb  = new TGeoBBox("ShADFlatBarb"  , ds    ,  7/2  , 8/2   ); 
  TGeoArb8 * shADFlatBara  = new TGeoArb8("shADFlatBara"  , 10/2. ); 
  TGeoArb8 * shADFlatBara2 = new TGeoArb8("shADFlatBara2" , 10/2. ); 
  TGeoTube * ShADBolts     = new TGeoTube("ShADBolts"     , 0.    ,  0.7  , 5.4/2 ); 
      for (Int_t i=0; i<4; i++)
  {
      shADFlatBara  -> SetVertex(i   , trFx[i]  , trFy[i]  ); 
      shADFlatBara  -> SetVertex(i+4 , trFx[i]  , trFy[i]  ); 
      shADFlatBara2 -> SetVertex(i   , trFx2[i] , trFy2[i] ); 
      shADFlatBara2 -> SetVertex(i+4 , trFx2[i] , trFy2[i] ); 
  }

  TGeoRotation * Rx45  = new TGeoRotation("Rx45" ,   0.,   0.,  45.) ;
  TGeoRotation * Rx45m = new TGeoRotation("Rx45m",   0.,   0., -45.) ;
  (new TGeoCombiTrans("trFlatBar1", -25.,  35/2., 0., Rx45m)) -> RegisterYourself();
  (new TGeoCombiTrans("trFlatBar2",  25., -35/2., 0., Rx45m)) -> RegisterYourself();
  (new TGeoCombiTrans("trFlatBar3",  25.,  35/2., 0,  Rx45 )) -> RegisterYourself();
  (new TGeoCombiTrans("trFlatBar4", -25., -35/2., 0., Rx45 )) -> RegisterYourself();
  TGeoCompositeShape * shFlatBarc = new TGeoCompositeShape("shFlatBarc", 
      "ShADFlatBar - ShADFlatBarb:trFlatBar4 - ShADFlatBarb:trFlatBar3 - ShADFlatBarb:trFlatBar2 - ShADFlatBarb:trFlatBar1 - shADFlatBara - shADFlatBara2");
  TGeoVolume * voADFlatBar = new TGeoVolume("voADFlatBar", shFlatBarc, kMedSteelSh);
  voADFlatBar -> SetLineColor(kGray);
  TGeoVolume * voADBolt   = new TGeoVolume("voADBolt", ShADBolts, kMedSteelSh);
  voADBolt -> SetLineColor(kGray);
  
  /////////////////////////////////////////////////////////////////////////////
  //
  // Make Vacuum Suport Vab 15 in RB24/1 FP (B-Side) /marvin.ascencio@pucp.edu.oe
  // Sheet (15x200x220)
  // LHC Drawing: LHCVHL__0024
  
  TGeoBBox * ShADSheet = new TGeoBBox("ShADSheet", 11, 10, 1.5/2);
  TGeoBBox * ShADSheetb = new TGeoBBox("ShADSheetb", ds, 7/2, 8/2); 

  (new TGeoCombiTrans("trSheet1", -11.,  10., 0., Rx45m)) -> RegisterYourself();
  (new TGeoCombiTrans("trSheet2",  11., -10., 0., Rx45m)) -> RegisterYourself();
  (new TGeoCombiTrans("trSheet3",  11.,  10., 0.,  Rx45 )) -> RegisterYourself();
  (new TGeoCombiTrans("trSheet4", -11., -10., 0., Rx45 )) -> RegisterYourself();
  TGeoCompositeShape * shSheet = new TGeoCompositeShape("shSheet", 
      "ShADSheet - ShADSheetb:trSheet4 - ShADSheetb:trSheet3 - ShADSheetb:trSheet2 - ShADSheetb:trSheet1");
  TGeoVolume * voADSheet = new TGeoVolume("voADSheet", shSheet, kMedSteelSh);
  voADSheet -> SetLineColor(kGray);
  //voWarmModuleSupport ->AddNode(voADSheet,2, new TGeoTranslation(0., 14.3,78.65));

  /////////////////////////////////////////////////////////////////////////////
  //
  // Make Vacuum Suport Vab 15 in RB24/1 FP (B-Side) /marvin.ascencio@pucp.edu.pe
  // TUBE epr.5 (100x80x700)
  // LHC Drawing: LHCVHL__0024
  // FOOT UPPER PLATE
  
  TGeoBBox * ShADTubeepr  = new TGeoBBox("ShADTubeepr"  , 5.0 , 4.  , 39.0 ); 
  TGeoBBox * ShADTubeeprb = new TGeoBBox("ShADTubeeprb" , 4.0 , 3.0 , 39.0     ); 
  TGeoCompositeShape * shTubeepr = new TGeoCompositeShape("shTubeepr", "ShADTubeepr - ShADTubeeprb");
  TGeoVolume * voADTubeepr = new TGeoVolume("voADTubeepr", shTubeepr, kMedSteelSh);
  voADTubeepr -> SetLineColor(kGreen);
  TGeoVolumeAssembly * voADsuppVab15 = new TGeoVolumeAssembly("voADsuppVab15");
  voADsuppVab15 -> AddNode(voADFlatBar,1 , new TGeoCombiTrans(14.3, 10.0+8.8, 0., Rz90 ));
  voADsuppVab15 -> AddNode(voADSheet, 1, new TGeoTranslation(14.3, 0., 86.65) );
  voADsuppVab15 -> AddNode(voADSheet, 2, new TGeoTranslation(14.3, 0., 79.75));
  voADsuppVab15 -> AddNode(voADTubeepr, 1, new TGeoTranslation(14.3, 0., 40.));
  voADsuppVab15 -> AddNode(voADBolt, 1, new TGeoTranslation(14.3+7.7, 8.3, 83.2));
  voADsuppVab15 -> AddNode(voADBolt, 2, new TGeoTranslation(14.3-7.7, 8.3, 83.2));
  voADsuppVab15 -> AddNode(voADBolt, 3, new TGeoTranslation(14.3+7.7, -8.3, 83.2));
  voADsuppVab15 -> AddNode(voADBolt, 4, new TGeoTranslation(14.3-7.7, -8.3, 83.2));

  ST0521057_VAMPA ->AddNode(RB241_WarmModule, 1, new TGeoCombiTrans(0., 0., 0.0, Rz90 ));
  ST0521057_VAMPA ->AddNode(voADsuppVab15, 1, new TGeoCombiTrans(-15.0+1, -119-1.5, -10.0 -1.5,Rx90m));
  return (TGeoVolume*) ST0521057_VAMPA;
}





