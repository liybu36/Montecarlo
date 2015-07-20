// --------------------------------------------------------------------------//
/** 
 * AUTHOR: Davide Franco
 * 
 */
// --------------------------------------------------------------------------//

#ifndef DSParametersH_H
#define DSParametersH_H
#include "G4ThreeVector.hh"
#include <vector>
#include "math.h"
#include <iostream>
#include "G4MaterialPropertiesTable.hh"
#include "G4EmSaturation.hh"

//---------------------------------------------------------------------------//

using namespace std;
class DSParameters  {
  private:
    DSParameters();
  public:
    static DSParameters* Get();
    
    virtual ~DSParameters() {}

    void       ReadDetectorGeometry();
    void       ReadDetectorProperties();
    void       ReadVetoPMTGeometry();
    void       ReadOpticalProperties();
    void       ReadOpticsTuning();
 
    G4double   GetTPCQE(G4double);
    G4double   GetWorldSizeX()                   {return fWorldSizeX;}        
    G4double   GetWorldSizeY()                   {return fWorldSizeY;}        
    G4double   GetWorldSizeZ()                   {return fWorldSizeZ;}        

   //Tank
    G4double   GetTankRmax()                     {return fTankRmax;}        
    G4double   GetTankHeight()                   {return fTankHeight;}	   
 
    G4double   GetSSSRadius()                    { return fSSSRadius; }

    // Cryostat
    G4double   GetCryostatShiftZ()               { return fCryostatShiftZ; }

    G4double   GetTrunkShiftZ()                  { return fTrunkShiftZ; }
    G4double   GetTrunkDiameter()                { return fTrunkDiameter; }
    G4double   GetTrunkThickness()               { return fTrunkThickness; }
    G4double   GetTrunkTopBottomOffsetZ()        { return fTrunkTopBottomOffsetZ; } 
    G4double   GetTrunkTopBottomOffsetX()        { return fTrunkTopBottomOffsetX; } 
    G4double   GetTrunkBottomHeight()            { return fTrunkBottomHeight; }
    G4double   GetTrunkMiddleHeight()            { return fTrunkMiddleHeight; }
    G4double   GetTrunkDistToCryostatAxis()      { return fTrunkDistToCryostatAxis; }
                                         
    // TPC                               
    G4double   GetTPCShiftZ()                    { return fTPCShiftZ; }
                                         
    G4double   GetTeflonSupportDiameter()        { return fTeflonSupportDiameter; }
    G4double   GetTeflonSupportThickness()       { return fTeflonSupportThickness; }
    G4double   GetTeflonSupportHeight()          { return fTeflonSupportHeight; }
                                        
    G4double   GetPMTAssemblyHeight()            { return fPMTAssemblyHeight; }
    G4double   GetTeflonCapDiameter()            { return fTeflonCapDiameter; }
    G4double   GetTeflonCapHeight()              { return fTeflonCapHeight; }
                                         
    G4double   GetLArBottomLayerThickness()      { return fLArBottomLayerThickness; }
    G4double   GetCathodeWindowDiameter()        { return fCathodeWindowDiameter; }
    G4double   GetCathodeWindowHeight()          { return fCathodeWindowHeight; }
    G4double   GetITOThickness()                 { return fITOThickness; }
                                         
    G4double   GetReflectorInnerDiameter()       { return fReflectorInnerDiameter; }
    G4double   GetReflectorOuterDiameter()       { return fReflectorOuterDiameter; }
    G4double   GetReflectorHeight()              { return fReflectorHeight; }
    G4double   GetAboveGridInnerDiameter()       { return fAboveGridInnerDiameter; }
    G4double   GetAboveGridOuterDiameter()       { return fAboveGridOuterDiameter; }
    G4double   GetAboveGridHeight()              { return fAboveGridHeight; }
    G4double   GetGasPocketThickness()           { return fGasPocketThickness; }
    G4double   GetTPBThickness()                 { return fTPBThickness;  }
                                         
    G4double   GetDivingBellHeight()             { return fDivingBellHeight; }
    G4double   GetDivingBellTopHeight()          { return fDivingBellTopHeight; }
    G4double   GetDivingBellOuterDiameter()      { return fDivingBellOuterDiameter; }
    G4double   GetDivingBellInnerDiameter()      { return fDivingBellInnerDiameter; }

    G4double   GetFieldRingsHeight()             { return fFieldRingsHeight; }
    G4double   GetFieldRingsThickness()          { return fFieldRingsThickness; }

    // TPC PMTs
    G4double   GetPMTBodyDiameter()              { return fPMTBodyDiameter; }
    G4double   GetPMTBodyHeight()                { return fPMTBodyHeight; }
    G4double   GetPMTHeadDiameter()              { return fPMTHeadDiameter; }
    G4double   GetPMTHeadHeight()                { return fPMTHeadHeight; } 
    G4double   GetPMTWallThickness()             { return fPMTWallThickness; }
    G4double   GetPMTWindowThickness()           { return fPMTWindowThickness; }
    G4double   GetPMTSpacing()                   { return fPMTSpacing; }
    G4double   GetPMTOffset()                    { return fPMTOffset; }

    // Veto PMT Positions
    G4int GetNVPMTs()                            { return fVPMTNum.size(); }
    vector<G4double>   GetVPMTNumber()           { return fVPMTNum; }
    vector<G4double>   GetVPMTTheta()            { return fVPMTTheta; }
    vector<G4double>   GetVPMTPhi()              { return fVPMTPhi; }

    // Properties
    G4double          GetPmtMaxQe()                     { return fPmtMaxQe;          }
    G4double          GetRealPmtMaxQe()                 { return fRealPmtMaxQe;      }
    G4double          GetPmtMaxQe_adjusted()            { return fPmtMaxQe_adjusted; }
    vector<G4double>  GetPmtQeWl()                      { return fPmtQeWl;           }
    vector<G4double>  GetPmtQe()                        { return fPmtQe;             }
    G4double          GetVPmtMaxQe()                    { return fVPmtMaxQe;         }
    G4double          GetVPmtMaxQe_adjusted()           { return fVPmtMaxQe_adjusted;}
    vector<G4double>  GetVPmtQeWl()                     { return fVPmtQeWl;          }
    vector<G4double>  GetVPmtQe()                       { return fVPmtQe;            }

    std::vector<G4double> ConvertTo4DoubleVector(const char* st);
    
    
    // Optics Tuning
    G4double  GetGridSteelRindScale()            { return fGridSteelRindScale;    }
    G4double  GetPhotocathodeUVRind()            { return fPhotocathodeUVRind; }
    G4double  GetPhotocathodeVisRind()           { return fPhotocathodeVisRind;     }
    G4double  GetFusedSilicaUVRind()             { return fFusedSilicaUVRind;	    }
    G4double  GetFusedSilicaVisRind()            { return fFusedSilicaVisRind;   }
    G4double  GetTPBUVRind()                     { return fTPBUVRind;}
    G4double  GetTPBVisRind()                    { return fTPBVisRind;  	}
    G4double  GetGaseousArgonUVAbs()             { return fGaseousArgonUVAbs;  	  }  
    G4double  GetGaseousArgonVisAbs()            { return fGaseousArgonVisAbs;    }
    G4double  GetLiquidArgonUVAbs()              { return fLiquidArgonUVAbs; }
    G4double  GetLiquidArgonVisAbs()             { return fLiquidArgonVisAbs; 	  }
    G4double  GetGridSteelUVAbs()                { return fGridSteelUVAbs;	  }
    G4double  GetGridSteelVisAbs()               { return fGridSteelVisAbs;       }
    G4double  GetTPBUVAbs()                      { return fTPBUVAbs ;}
    G4double  GetTPBVisAbs()                     { return fTPBVisAbs;  	}
    G4double  GetFusedSilicaUVAbs()              { return fFusedSilicaUVAbs;	  }
    G4double  GetFusedSilicaVisAbs()             { return fFusedSilicaVisAbs;         }
    G4double  GetWLSAbsorptionFactor()           { return fWLSAbsorptionFactor ;}
    G4double  GetWLSMeanNumberPhotons()          { return fWLSMeanNumberPhotons;   }
    G4double  GetWLSTimeConstant_ns()            { return fWLSTimeConstant_ns;       }
    G4double  GetWLSEfficiency()                 { return fWLSEfficiency;	  }
    G4int     GetWithITO()                       { return fWithITO;          }
    G4int     GetWithGasPocket()                 { return fWithGasPocket; }
    G4int     GetWithNewGridModel()              { return fWithNewGridModel; }
    G4double  GetGridNormalTransparency()        { return fGridNormalTransparency;	}
    G4double  GetGArRindexScale()                { return fGArRindexScale; }
    G4double  GetLArRayleighScale()              { return fLArRayleighScale; 	  }
    G4double  GetFSilicaRaylVisLength()          { return fFSilicaRaylVisLength; }
    G4double  GetFSilicaRaylUVLength()           { return fFSilicaRaylUVLength; }
    G4double  GetLArGridUVRef()                  { return fLArGridUVRef;	}
    G4double  GetLArGridVisRef()                 { return fLArGridVisRef;	}
    G4double  GetTeflonTPBUVRef()                { return fTeflonTPBUVRef;	}
    G4double  GetTeflonTPBVisRef()               { return fTeflonTPBVisRef;	}
    G4double  GetPMTLArUVRef()                   { return fPMTLArUVRef;		  }
    G4double  GetPMTLArVisRef()                  { return fPMTLArVisRef;	}
    G4double  GetTeflonLArUVRef()                { return fTeflonLArUVRef;	}
    G4double  GetTeflonLArVisRef()               { return fTeflonLArVisRef;	}
    G4double  GetArTPBVisTran()                  { return fArTPBVisTran;	}
    G4double  GetPArRind()                       { return fPArRind;        }
    
    
    
    /*
    // Optical Properties
    G4MaterialPropertiesTable *GetLumirrorMPT() { return fLumirrorMPT; }
    G4MaterialPropertiesTable *GetElectropolishedStainlessSteelMPT() { return fElectropolishedStainlessSteelMPT; }
    G4MaterialPropertiesTable *GetUntreatedStainlessSteelMPT() { return fUntreatedStainlessSteelMPT; }
    G4MaterialPropertiesTable *GetAluminumFoilMPT() { return fAluminumFoilMPT; }
    G4MaterialPropertiesTable *GetPMTBackMPT() { return fPMTBackMPT; }
    G4MaterialPropertiesTable *GetVPhotocathodeMPT() { return fVPhotocathodeMPT; }
  */

    // Optical Properties
    G4int     GetLumirrorNumEntries() { return fLumirror_Energy.size(); }
    G4int     GetVQENumEntries() { return  fVPhotocathode_Energy.size(); }
    G4int     GetQENumEntries() { return fPhotocathode_Energy.size();}
    vector<G4double> GetLumirrorEnergy() { return fLumirror_Energy; }
    vector<G4double> GetLumirrorSpecularLobe() { return fLumirror_SpecularLobe; }
    vector<G4double> GetLumirrorSpecularSpike() { return fLumirror_SpecularSpike; }
    vector<G4double> GetLumirrorBackscatter() { return fLumirror_Backscatter; }
    vector<G4double> GetLumirrorReflectance() { return fLumirror_Reflectance; }
    vector<G4double> GetLumirrorEfficiency() { return fLumirror_Efficiency; }
    vector<G4double> GetLumirrorRindex() { return fLumirror_Rindex; }
    vector<G4double> GetElectropolishedSSEnergy() { return fElectropolishedStainlessSteel_Energy; }
    vector<G4double> GetElectropolishedSSSpecularLobe() { return fElectropolishedStainlessSteel_SpecularLobe; }
    vector<G4double> GetElectropolishedSSSpecularSpike() { return fElectropolishedStainlessSteel_SpecularSpike; }
    vector<G4double> GetElectropolishedSSBackscatter() { return fElectropolishedStainlessSteel_Backscatter; }
    vector<G4double> GetElectropolishedSSReflectance() { return fElectropolishedStainlessSteel_Reflectance; }
    vector<G4double> GetElectropolishedSSEfficiency() { return fElectropolishedStainlessSteel_Efficiency; }
    vector<G4double> GetElectropolishedSSRindex() { return fElectropolishedStainlessSteel_Rindex; }
    vector<G4double> GetUntreatedSSEnergy() { return fUntreatedStainlessSteel_Energy; }
    vector<G4double> GetUntreatedSSSpecularLobe() { return fUntreatedStainlessSteel_SpecularLobe; }
    vector<G4double> GetUntreatedSSSpecularSpike() { return fUntreatedStainlessSteel_SpecularSpike; }
    vector<G4double> GetUntreatedSSBackscatter() { return fUntreatedStainlessSteel_Backscatter; }
    vector<G4double> GetUntreatedSSReflectance() { return fUntreatedStainlessSteel_Reflectance; }
    vector<G4double> GetUntreatedSSEfficiency() { return fUntreatedStainlessSteel_Efficiency; }
    vector<G4double> GetUntreatedSSRindex() { return fUntreatedStainlessSteel_Rindex; }
    vector<G4double> GetAluminumFoilEnergy() { return fAluminumFoil_Energy; }
    vector<G4double> GetAluminumFoilSpecularLobe() { return fAluminumFoil_SpecularLobe; }
    vector<G4double> GetAluminumFoilSpecularSpike() { return fAluminumFoil_SpecularSpike; }
    vector<G4double> GetAluminumFoilBackscatter() { return fAluminumFoil_Backscatter; }
    vector<G4double> GetAluminumFoilReflectance() { return fAluminumFoil_Reflectance; }
    vector<G4double> GetAluminumFoilEfficiency() { return fAluminumFoil_Efficiency; }
    vector<G4double> GetAluminumFoilRindex() { return fAluminumFoil_Rindex; }
    vector<G4double> GetVPhotocathodeEnergy() { return fVPhotocathode_Energy; }
    vector<G4double> GetVPhotocathodeSpecularLobe() { return fVPhotocathode_SpecularLobe; }
    vector<G4double> GetVPhotocathodeSpecularSpike() { return fVPhotocathode_SpecularSpike; }
    vector<G4double> GetVPhotocathodeBackscatter() { return fVPhotocathode_Backscatter; }
    vector<G4double> GetVPhotocathodeReflectance() { return fVPhotocathode_Reflectance; }
    vector<G4double> GetVPhotocathodeEfficiency() { return fVPhotocathode_Efficiency; }
    vector<G4double> GetVPhotocathodeRindex() { return fVPhotocathode_Rindex; }
    vector<G4double> GetPMTBackEnergy() { return fPMTBack_Energy; }
    vector<G4double> GetPMTBackSpecularLobe() { return fPMTBack_SpecularLobe; }
    vector<G4double> GetPMTBackSpecularSpike() { return fPMTBack_SpecularSpike; }
    vector<G4double> GetPMTBackBackscatter() { return fPMTBack_Backscatter; }
    vector<G4double> GetPMTBackReflectance() { return fPMTBack_Reflectance; }
    vector<G4double> GetPMTBackEfficiency() { return fPMTBack_Efficiency; }
    vector<G4double> GetPMTBackRindex() { return fPMTBack_Rindex; }
    vector<G4double> GetPhotocathodeEnergy() { return fPhotocathode_Energy; }
    vector<G4double> GetPhotocathodeSpecularLobe() { return fPhotocathode_SpecularLobe; }
    vector<G4double> GetPhotocathodeSpecularSpike() { return fPhotocathode_SpecularSpike; }
    vector<G4double> GetPhotocathodeBackscatter() { return fPhotocathode_Backscatter; }
    vector<G4double> GetPhotocathodeReflectance() { return fPhotocathode_Reflectance; }
    vector<G4double> GetPhotocathodeEfficiency() { return fPhotocathode_Efficiency; }
    vector<G4double> GetPhotocathodeRindex() { return fPhotocathode_Rindex; }

    // Cryostat Sheath Properties
    G4double*                   GetCryoSheathR() { return fCryoSheathR; }
    G4double*                   GetCryoSheathZ() { return fCryoSheathZ; }

    // EM Saturation
    void AddSaturation(G4EmSaturation* sat) { emSaturation = sat; }
    void RemoveSaturation() { emSaturation = NULL; }
    G4EmSaturation *GetSaturation() { return emSaturation; }

  private:
  
    static DSParameters *me;       
    
    G4String   fDSGeometryFileName ;
    
    float      fWorldSizeX;
    float      fWorldSizeY;
    float      fWorldSizeZ;
    float      fTankRmax;
    float      fTankHeight;
    
    G4double   fSSSRadius;

    // Cryostat
    G4double   fCryostatShiftZ;

    G4double   fTrunkShiftZ;
    G4double   fTrunkDiameter;
    G4double   fTrunkThickness;
    G4double   fTrunkTopBottomOffsetZ; 
    G4double   fTrunkTopBottomOffsetX; 
    G4double   fTrunkBottomHeight;
    G4double   fTrunkMiddleHeight;
    G4double   fTrunkDistToCryostatAxis;

    // TPC
    G4double   fTPCShiftZ;

    G4double   fTeflonSupportDiameter;
    G4double   fTeflonSupportThickness;
    G4double   fTeflonSupportHeight;

    G4double   fPMTAssemblyHeight;
    G4double   fTeflonCapDiameter;
    G4double   fTeflonCapHeight;

    G4double   fLArBottomLayerThickness;
    G4double   fCathodeWindowDiameter;
    G4double   fCathodeWindowHeight;
    G4double   fITOThickness;

    G4double   fReflectorInnerDiameter;
    G4double   fReflectorOuterDiameter;
    G4double   fReflectorHeight;
    G4double   fAboveGridInnerDiameter;
    G4double   fAboveGridOuterDiameter;
    G4double   fAboveGridHeight;
    G4double   fGasPocketThickness;
    G4double   fTPBThickness;

    G4double   fDivingBellHeight;
    G4double   fDivingBellTopHeight;
    G4double   fDivingBellOuterDiameter;
    G4double   fDivingBellInnerDiameter;

    G4double   fFieldRingsHeight;
    G4double   fFieldRingsThickness;

    // TPC PMTs
    G4double   fPMTBodyDiameter;
    G4double   fPMTBodyHeight;
    G4double   fPMTHeadDiameter;
    G4double   fPMTHeadHeight;
    G4double   fPMTWallThickness;
    G4double   fPMTWindowThickness;
    G4double   fPMTSpacing;
    G4double   fPMTOffset;

    // Veto PMT Positions
    //static const int fNVPMTs = 110;
    vector<G4double> fVPMTNum;
    vector<G4double> fVPMTTheta;
    vector<G4double> fVPMTPhi;
    
    
   // Optics Tuning
    G4double   fGridSteelRindScale;
    G4double   fPhotocathodeUVRind;
    G4double   fPhotocathodeVisRind;
    G4double   fFusedSilicaUVRind;
    G4double   fFusedSilicaVisRind;
    G4double   fTPBUVRind;
    G4double   fTPBVisRind;
    G4double   fGaseousArgonUVAbs;
    G4double   fGaseousArgonVisAbs;
    G4double   fLiquidArgonUVAbs;
    G4double   fLiquidArgonVisAbs;
    G4double   fGridSteelUVAbs ;
    G4double   fGridSteelVisAbs;
    G4double   fFusedSilicaUVAbs;
    G4double   fFusedSilicaVisAbs;
    G4double   fTPBUVAbs;
    G4double   fTPBVisAbs;
    G4double   fWLSAbsorptionFactor;       
    G4double   fWLSMeanNumberPhotons;      
    G4double   fWLSTimeConstant_ns;        
    G4double   fWLSEfficiency;             
    G4int      fWithITO;                   
    G4int      fWithGasPocket;             
    G4int      fWithNewGridModel;  
    G4double   fGridNormalTransparency;           
    G4double   fGArRindexScale;          
    G4double   fLArRayleighScale; 
    G4double   fFSilicaRaylVisLength;
    G4double   fFSilicaRaylUVLength;
    G4double   fLArGridUVRef;              
    G4double   fLArGridVisRef;             
    G4double   fTeflonTPBUVRef;            
    G4double   fTeflonTPBVisRef;           
    G4double   fPMTLArUVRef;               
    G4double   fPMTLArVisRef;              
    G4double   fTeflonLArUVRef;            
    G4double   fTeflonLArVisRef;           
    G4double   fArTPBVisTran;           
    G4double   fPArRind;

    // Optical Boundary Properties
    
    G4int fLumirrorSkipNEntries ;    // initialized in the DSParameters constructor 
    //static const G4int fLumirrorNumEntries = (int)601/fLumirrorSkipNEntries;
    //const G4int fVQENumEntries = 40;
    //const G4int fQENumEntries = 59;
    // Energy vectors
    vector<G4double> fLumirror_Energy;
    vector<G4double> fElectropolishedStainlessSteel_Energy;
    vector<G4double> fUntreatedStainlessSteel_Energy;
    vector<G4double> fAluminumFoil_Energy;
    vector<G4double> fPMTBack_Energy;
    vector<G4double> fVPhotocathode_Energy;
    vector<G4double> fPhotocathode_Energy;
    // Reflectance vectors
    vector<G4double> fLumirror_Reflectance;
    vector<G4double> fElectropolishedStainlessSteel_Reflectance;
    vector<G4double> fUntreatedStainlessSteel_Reflectance;
    vector<G4double> fAluminumFoil_Reflectance;
    vector<G4double> fPMTBack_Reflectance;
    vector<G4double> fVPhotocathode_Reflectance;
    vector<G4double> fPhotocathode_Reflectance;
    // Efficiency vectors
    vector<G4double> fLumirror_Efficiency;
    vector<G4double> fElectropolishedStainlessSteel_Efficiency;
    vector<G4double> fUntreatedStainlessSteel_Efficiency;
    vector<G4double> fAluminumFoil_Efficiency;
    vector<G4double> fPMTBack_Efficiency;
    vector<G4double> fVPhotocathode_Efficiency;
    vector<G4double> fPhotocathode_Efficiency;
    // Specular lobe vectors
    vector<G4double> fLumirror_SpecularLobe;
    vector<G4double> fElectropolishedStainlessSteel_SpecularLobe;
    vector<G4double> fUntreatedStainlessSteel_SpecularLobe;
    vector<G4double> fAluminumFoil_SpecularLobe;
    vector<G4double> fPMTBack_SpecularLobe;
    vector<G4double> fVPhotocathode_SpecularLobe;
    vector<G4double> fPhotocathode_SpecularLobe;
    // Specular spike vectors
    vector<G4double> fLumirror_SpecularSpike;
    vector<G4double> fElectropolishedStainlessSteel_SpecularSpike;
    vector<G4double> fUntreatedStainlessSteel_SpecularSpike;
    vector<G4double> fAluminumFoil_SpecularSpike;
    vector<G4double> fPMTBack_SpecularSpike;
    vector<G4double> fVPhotocathode_SpecularSpike;
    vector<G4double> fPhotocathode_SpecularSpike;
    // Backscatter vectors
    vector<G4double> fLumirror_Backscatter;
    vector<G4double> fElectropolishedStainlessSteel_Backscatter;
    vector<G4double> fUntreatedStainlessSteel_Backscatter;
    vector<G4double> fAluminumFoil_Backscatter;
    vector<G4double> fPMTBack_Backscatter;
    vector<G4double> fVPhotocathode_Backscatter;
    vector<G4double> fPhotocathode_Backscatter;
    // Rindex vectors
    vector<G4double> fLumirror_Rindex;
    vector<G4double> fElectropolishedStainlessSteel_Rindex;
    vector<G4double> fUntreatedStainlessSteel_Rindex;
    vector<G4double> fAluminumFoil_Rindex;
    vector<G4double> fPMTBack_Rindex;
    vector<G4double> fVPhotocathode_Rindex;
    vector<G4double> fPhotocathode_Rindex;
  /*
    // Material properties table
    G4MaterialPropertiesTable *fLumirrorMPT;
    G4MaterialPropertiesTable *fElectropolishedStainlessSteelMPT;
    G4MaterialPropertiesTable *fUntreatedStainlessSteelMPT;
    G4MaterialPropertiesTable *fAluminumFoilMPT;
    G4MaterialPropertiesTable *fVPhotocathodeMPT;
    G4MaterialPropertiesTable *fPMTBackMPT;
  */
    // Properties
    G4double   fPmtMaxQe;
    G4double   fRealPmtMaxQe;
    G4double   fPmtMaxQe_adjusted;
    G4double   fVPmtMaxQe;
    G4double   fVPmtMaxQe_adjusted;

    vector<G4double>   fPmtQeWl;
    vector<G4double>   fPmtQe;
    vector<G4double>   fVPmtQeWl;
    vector<G4double>   fVPmtQe;

    // Dimensions for the lumirror sheath around the outside of the cryostat
    G4double fCryoSheathR[2];
    G4double fCryoSheathZ[2];

    G4EmSaturation* emSaturation;
};



#endif
/*
 * $Log: DSParameters.hh,v $
 * Revision 1.7  2015/01/17 11:31:45  pagnes
 * PAr model added form optical tuning
 *
 * Revision 1.6  2015/01/07 16:45:59  pagnes
 * changed veto optical properties format from arrays to vectors
 *
 * Revision 1.5  2014/11/05 15:47:09  pagnes
 * temporary optics tuning
 *
 * Revision 1.4  2014/10/13 18:43:46  swesterd
 * fixed veto PMT positions
 *
 * Revision 1.3  2014/07/16 08:23:12  pagnes
 * QE scaling to 1.0 added (/ds/manager/fast_simulation xxx)
 *
 * Revision 1.2  2014/06/03 13:31:33  meyers
 * Migrate TPC grid and ITO optics updates to g4ds10
 *
 * Revision 1.1  2014/05/07 12:20:54  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.17  2014/03/19 16:37:36  dfranco
 * update external configuration reader for optics tuning
 *
 * Revision 1.16  2014/03/11 16:50:00  dfranco
 * read optical parameters from external file (DSOptics.dat)
 *
 * Revision 1.15  2013/08/27 07:13:12  swesterd
 * add visible energy for the neutron veto
 *
 * Revision 1.14  2013/06/19 18:35:25  swesterd
 * added DSScintCelll and made tpc PMTs' QE and reflections work like veto PMTs
 *
 * Revision 1.13  2013/06/05 23:03:28  swesterd
 * moved optical boundary MPTs to DSMaterial and gave the trunks optical boundary properties consistent with untreated stainless steel
 *
 * Revision 1.12  2013/06/04 23:38:46  swesterd
 * Changed the length of the Lumirror reflectance array to be an adjustable size, set it to 60 elements
 *
 * Revision 1.11  2013/06/04 14:11:38  dfranco
 * Added a function returning the TPC QE in DSParameters, applied in the tracking action
 *
 * Revision 1.10  2013/05/27 23:59:00  swesterd
 * added a (currently commented out) Lumirror sheath to the cryostat and introduced DSOpBoundaryProcess to try to figure out why the boundaries are being screwy, with some edits so that it can handle constant and vector properties with freaking out
 *
 * Revision 1.9  2013/05/25 07:58:22  swesterd
 * Got the veto PMT optical boundaries all working along with photocathode optical properties, added PMT quantum efficiency to DSTrackingAction, and added a function to DSTrackingAction that locates and quadratically interpolates points in data, for getting useful QEs
 *
 * Revision 1.8  2013/05/07 23:06:26  swesterd
 * added optical boundaries and Lumirror in the veto
 *
 * Revision 1.7  2013/05/06 14:59:52  perassos
 * Updates on the TPC surface properties and geometry
 *
 * Revision 1.6  2013/05/01 08:20:23  swesterd
 * added boron-loaded scintillator optical properties
 *
 * Revision 1.5  2013/04/26 16:17:19  perassos
 * TPC PMT QE added
 *
 * Revision 1.4  2013/04/19 16:10:09  perassos
 * Added ITO and TPB layers
 *
 * Revision 1.3  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
