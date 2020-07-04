/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.0 $
// $Date: 2019-08-22 $
                                                      
// Written: Kuanshi Zhong
// Created: 08/2019
//



#ifndef DuctileFracture_h
#define DuctileFracture_h

#include <UniaxialMaterial.h>

class DuctileFracture : public UniaxialMaterial
{
 public:

	 DuctileFracture(int tag, UniaxialMaterial &material,double Dmax = 1.0,
		 double c_mono = 0.20,double c_cycl = 0.15,double c_symm = 1.30,
		 double E_s = 29000.0,double minStrain = -1.0e16,double maxStrain =  1.0e16,
		 double esu = 1.0e16,double k1 = 1,double k2 = 0,double db = 0,double b1 = 1000,double b2 = -1000,
		 double c_dete = 0.0);
	 DuctileFracture();
	 ~DuctileFracture();

	  const char *getClassType(void) const {return "DuctileFracture";};
  
	  int setTrialStrain(double strain, double strainRate = 0.0); 
	  double getStrain(void);          
	  double getStrainRate(void);
	  double getStress(void);
	  double getTangent(void);
	  double getDampTangent(void);
	  double getInitialTangent(void) {return theMaterial->getInitialTangent();}
  
	  int commitState(void);
	  int revertToLastCommit(void);    
	  int revertToStart(void);        
  
	  UniaxialMaterial *getCopy(void);
  
	  int sendSelf(int commitTag, Channel &theChannel);  
	  int recvSelf(int commitTag, Channel &theChannel, 
			   FEM_ObjectBroker &theBroker);    
  
	  void Print(OPS_Stream &s, int flag =0);

	  Response *setResponse (const char **argv, int argc, OPS_Stream &s);
	  int getResponse (int responseID, Information &matInformation);    
	  bool hasFailed(void);  

 protected:
  
 private:
	  UniaxialMaterial *theMaterial;
  
	  double DI; // Damage index
	  double DI_VGM; // Void growth damage component
	  double DI_MVC; // Multi-void coalescence damage compoent
	  double ep_prev; // Previous plastic strain
	  double ep_curr; // Current plastic strain
	  double dep; // Incremental plastic strain
	  double cep_comp; // Cumulative compressive plastic strain
	  double es_local; // Local strain
	  double T; // Triaxiality
	  double es_max; // The maximum steel strain
	  double es_min; // The minimum steel strain
	  double e_memo; // The strain memory factor
  
	  double Dmax;
	  double c_mono;
	  double c_cycl;
	  double c_symm;
	  double E_s;
	  double minStrain;
	  double maxStrain;
	  double esu; // necking strain (esu = 0: no necking)
	  double k1; // necking local strain amplification (k1 = 1: no amplification) 
	  double k2; // necking local Triaxiality amplification (k2 = 0: no amplification)
	  double db; // bar-diameter for buckle amplification
	  double b1; // buckle initial-imperfection coefficient (b1 = 1000: no effect)
	  double b2; // buckle curvature coefficient (b2 = -1000: no local buckle)
	  double c_dete; // strength deterioration coefficient (c_dete = 0: no deterioration)
	  bool Cfailed;
	  double trialStrain;
  
};


#endif

