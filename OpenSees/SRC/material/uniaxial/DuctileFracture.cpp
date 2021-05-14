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
// Created: 07/2019
//

#include <math.h>

#include <stdlib.h>
#include <MaterialResponse.h>
#include <Information.h>
#include <string.h>

#include <DuctileFracture.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <OPS_Globals.h>
#include <OPS_Stream.h>
#include <elementAPI.h>

void* OPS_DuctileFracture()
{
	int numdata = OPS_GetNumRemainingInputArgs();
	if (numdata < 2) {
		opserr << "WARNING insufficient arguments\n";
		opserr << "Want: uniaxialMaterial DuctileFracture tag? matTag?";
		opserr << " <-D_max dmax?> <-c_mono c_mono?> <-c_cycl c_cycl?> <-c_symm c_symm?>" << endln;
		opserr << " <-E_s E_s> <-min min?> <-max max?>" << endln;
		opserr << " <-esu esu> <-k1 k1> <-k2 k2> <-db db> <-b1 b1> <-b2 b2> <-c_dete c_dete>" << endln;
		return 0;
	}

	int idata[2];
	numdata = 2;
	if (OPS_GetIntInput(&numdata, idata) < 0) {
		opserr << "WARNING invlid int inputs\n";
		return 0;
	}

	double Dmax = 1.0;
	double c_mono = 0.20;
	double c_cycl = 0.15;
	double c_symm = 1.30;
	double E_s = 29000;
	double epsmin = NEG_INF_STRAIN;
	double epsmax = POS_INF_STRAIN;
	double esu = POS_INF_STRAIN; // necking strain (esu = 0: no necking)
	double k1 = 1; // necking local strain amplification (k1 = 1: no amplification) 
	double k2 = 0; // necking local Triaxiality amplification (k2 = 0: no amplification)
	double db = 0; // bar-diameter for buckle amplification
	double b1 = 1000; // buckle initial-imperfection coefficient (b1 = 1000: no effect)
	double b2 = -1000; // buckle curvature coefficient (b2 = -1000: no local buckle)
	double c_dete = 0.0; // deteriroation coefficient (c_dete = 0: no deterioration)
	numdata = 1;

	while (OPS_GetNumRemainingInputArgs() > 1) {
		const char* type = OPS_GetString();
		if (strcmp(type, "-Dmax") == 0) {
			if (OPS_GetDouble(&numdata, &Dmax) < 0) {
				opserr << "WARNING invalid double inputs\n";
				return 0;
			}
		}
		else if (strcmp(type, "-c_mono") == 0) {
			if (OPS_GetDouble(&numdata, &c_mono) < 0) {
				opserr << "WARNING invalid double inputs\n";
				return 0;
			}
		}
		else if (strcmp(type, "-c_cycl") == 0) {
			if (OPS_GetDouble(&numdata, &c_cycl) < 0) {
				opserr << "WARNING invalid double inputs\n";
				return 0;
			}
		}
		else if (strcmp(type, "-c_symm") == 0) {
			if (OPS_GetDouble(&numdata, &c_symm) < 0) {
				opserr << "WARNING invalid double inputs\n";
				return 0;
			}
		}
		else if (strcmp(type, "-E_s") == 0) {
			if (OPS_GetDouble(&numdata, &E_s) < 0) {
				opserr << "WARNING invalid double inputs\n";
				return 0;
			}
		}
		else if (strcmp(type, "-min") == 0) {
			if (OPS_GetDouble(&numdata, &epsmin) < 0) {
				opserr << "WARNING invalid double inputs\n";
				return 0;
			}
		}
		else if (strcmp(type, "-max") == 0) {
			if (OPS_GetDouble(&numdata, &epsmax) < 0) {
				opserr << "WARNING invalid double inputs\n";
				return 0;
			}
		}
		else if (strcmp(type, "-esu") == 0) {
			if (OPS_GetDouble(&numdata, &esu) < 0) {
				opserr << "WARNING invalid double inputs\n";
				return 0;
			}
		}
		else if (strcmp(type, "-k1") == 0) {
			if (OPS_GetDouble(&numdata, &k1) < 0) {
				opserr << "WARNING invalid double inputs\n";
				return 0;
			}
		}
		else if (strcmp(type, "-k2") == 0) {
			if (OPS_GetDouble(&numdata, &k2) < 0) {
				opserr << "WARNING invalid double inputs\n";
				return 0;
			}
		}
		else if (strcmp(type, "-db") == 0) {
			if (OPS_GetDouble(&numdata, &db) < 0) {
				opserr << "WARNING invalid double inputs\n";
				return 0;
			}
		}
		else if (strcmp(type, "-b1") == 0) {
			if (OPS_GetDouble(&numdata, &b1) < 0) {
				opserr << "WARNING invalid double inputs\n";
				return 0;
			}
		}
		else if (strcmp(type, "-b2") == 0) {
			if (OPS_GetDouble(&numdata, &b2) < 0) {
				opserr << "WARNING invalid double inputs\n";
				return 0;
			}
		}
		else if (strcmp(type, "-c_dete") == 0) {
			if (OPS_GetDouble(&numdata, &c_dete) < 0) {
				opserr << "WARNING invalid double inputs\n";
				return 0;
			}
		}
		

		UniaxialMaterial* mat = OPS_getUniaxialMaterial(idata[1]);
		if (mat == 0) {
			opserr << "WARNING component material does not exist\n";
			opserr << "Component material: " << idata[1];
			opserr << "\nuniaxialMaterial DuctileFracture: " << idata[0] << endln;
			return 0;
		}

		UniaxialMaterial* theMat = new DuctileFracture(idata[0], *mat,
			Dmax, c_mono, c_cycl, c_symm, E_s, epsmin, epsmax, 
			esu, k1, k2, db, b1, b2,c_dete);
		if (theMat == 0) {
			opserr << "WARNING: failed to create DuctileFracture material\n";
			return 0;
		}

		return theMat;
	}
}

DuctileFracture::DuctileFracture(int tag,UniaxialMaterial &material,
	double dmax, double epf, double lambda, double beta,
	double Es, double epsmin, double epsmax, double e_su, double k_1, double k_2,
	double d_b, double b_1, double b_2,double cdete)
	:UniaxialMaterial(tag,MAT_TAG_DuctileFracture), theMaterial(0), 
	Cfailed(false), trialStrain(0)
{
  DI  = 0; //Damage index
  DI_VGM	= 0; // Void growth damage component
  DI_MVC	= 0; // Multi-void coalescence damage compoent
  ep_prev	= 0; // Previous plastic strain
  ep_curr	= 0; // Current plastic strain
  dep		= 0; // Incremental plastic strain
  cep_comp	= 0; // Cumulative compressive plastic strain
  es_local = 0; // Local strain
  T = 0; // Triaxiality
  es_max = 0; // The maximum steel strain
  es_min = 0; // The minimum steel strain
  e_memo = 0; // The strain memory factor

  if ( dmax > 10.0 || dmax < 0.0 ) {
    opserr << "DuctileFracture::DuctileFracture " <<
      "- Dmax must be between 0 and 10, assuming Dmax = 1\n" ;
    Dmax    = 1.0;
  } else 
    Dmax    = dmax;
  
  c_mono	= epf;
  c_cycl	= lambda;
  c_symm	= beta;
  E_s		= Es;
  minStrain = epsmin;
  maxStrain = epsmax;
  esu = e_su;
  k1 = k_1;
  k2 = k_2;
  db = d_b;
  b1 = b_1;
  b2 = b_2;
  c_dete = cdete;
  
  theMaterial = material.getCopy();
  if (theMaterial == 0) {
    opserr <<  "DuctileFracture::DuctileFracture " <<
      " -- failed to get copy of material\n" ;
    exit(-1);
  }
}

DuctileFracture::DuctileFracture()
	:UniaxialMaterial(0, MAT_TAG_DuctileFracture), theMaterial(0),
	Cfailed(false), trialStrain(0)
{
	DI = 0; //Damage index
	DI_VGM = 0; // Void growth damage component
	DI_MVC = 0; // Multi-void coalescence damage compoent
	ep_prev = 0; // Previous plastic strain
	ep_curr = 0; // Current plastic strain
	dep = 0; // Incremental plastic strain
	cep_comp = 0; // Cumulative compressive plastic strain
	es_local = 0; // Local strain
	T = 0; // Triaxiality
	es_max = 0; // The maximum steel strain
	es_min = 0; // The minimum steel strain
	e_memo = 0; // The strain memory factor

	Dmax = 1.0;
	c_mono = 0;
	c_cycl = 0;
	c_symm = 0;
	E_s = 0;
	minStrain = NEG_INF_STRAIN;
	maxStrain = POS_INF_STRAIN;
	esu = 0;
	k1 = 0;
	k2 = 0;
	db = 0;
	b1 = 0;
	b2 = 0;
	c_dete = 0;
}

DuctileFracture::~DuctileFracture()
{
  if (theMaterial)
    delete theMaterial;
}

int 
DuctileFracture::setTrialStrain(double strain, double strainRate)
{
  if (Cfailed) {
    trialStrain = strain;
    // return 0;
    return theMaterial->setTrialStrain(strain, strainRate);
  } else {
    Cfailed = false;
    trialStrain = strain;
    return theMaterial->setTrialStrain(strain, strainRate);
  }
}

double 
DuctileFracture::getStress(void)
{
  double modifier = 1.0;
  double damageloc = 1.0-Dmax+DI;
  if (Cfailed)
	  // Reduce stress to 0.0 
	  return theMaterial->getStress()*1.0e-8;
  else if (DI_MVC > 1.0) {
	  modifier = 1.0 / sqrt(pow(DI_MVC, c_dete));
	  return theMaterial->getStress()*modifier;
  }
  else
    return theMaterial -> getStress();
}

double 
DuctileFracture::getTangent(void)
{
  double modifier = 1.0;
  double damageloc = 1.0-Dmax+DI;
  if (Cfailed)
    // Reduce tangent to 0.0 
    return 1.0e-8*theMaterial->getInitialTangent();
  else if (DI_MVC > 1.0) {
	  modifier = 1.0 / sqrt(pow(DI_MVC, c_dete));
	  return theMaterial->getTangent()*modifier;
  }
  else
    return theMaterial->getTangent()*modifier;
}

double 
DuctileFracture::getDampTangent(void)
{
  if (Cfailed)
    return 0.0;
  else
    return theMaterial->getDampTangent();
}

double 
DuctileFracture::getStrain(void)
{
  return theMaterial->getStrain();
}

double 
DuctileFracture::getStrainRate(void)
{
  return theMaterial->getStrainRate();
}

int 
DuctileFracture::commitState(void)
{	
  if (Cfailed) {
    return 0;
  }

  if (trialStrain >= maxStrain || trialStrain <= minStrain) { 
      Cfailed = true;
      DI = Dmax;
      return 0;
  }
  
  // Update es_max
  if (trialStrain > es_max) {
	  es_max = trialStrain;
  }

  // Update es_min
  if (trialStrain < es_min) {
	  es_min = trialStrain;
  }

  // strain memory for c_symm and c_cycl
  e_memo = fabs(es_max - es_min) / 0.05;
  if (e_memo > 1.0) {
	  e_memo = 1.0;
  }

  // Get current stress
  double s_curr = theMaterial->getStress();

  // Necking amplification model
  if (trialStrain > esu) {
	  // Apply necking coefficient
	  es_local = esu + k1 * (trialStrain - esu);
	  T = 0.33 + k2 * (trialStrain - esu);
  }
  else {
	  es_local = trialStrain;
	  T = 0.33;
  }

  // Buckle amplification model
  es_local = es_local - 0.5 * db * (b1 * sinh((es_max - trialStrain) / b2));

  // Compute current plastic strain
  ep_curr = es_local-s_curr/E_s;
  
  // Compute current incremental plastic strain
  dep = ep_curr - ep_prev;  
  
  // Constant triaxiality T = 0.33
 
  // Compute damage index components
  if (dep > 0) {
	  DI_VGM = DI_VGM+c_mono*(((c_symm-1.0)*e_memo+1.0)*exp(1.3*T)-exp(-1.3*T))*fabs(dep);
  } else if (dep < 0) {
	  DI_VGM = DI_VGM+c_mono*(((c_symm-1.0)*e_memo+1.0)*exp(-1.3*T)-exp(1.3*T))*fabs(dep);
	  if (DI_VGM < 0) {
		  DI_VGM = 0;
	  }
	  cep_comp = cep_comp+fabs(dep);
  }
  // Move the multi-void coalescence outside (KZ - 05/13/21)
  DI_MVC = exp(c_cycl*e_memo*cep_comp);
  
  // Compute damage index
  DI = DI_VGM*DI_MVC;

  // Flag failure if we have reached that point
  if (DI >= Dmax ) {
	  Cfailed = true;
	  opserr << "DuctileFracture: material tag " << this->getTag() << " failed\n";
	  } else {
      Cfailed = false;
    }
  
  // Update previous plastic strain
  ep_prev = ep_curr;
  
  // Check if failed at current step
  if (Cfailed) {
    return 0;
  }
  else 
    return theMaterial->commitState();
}

int 
DuctileFracture::revertToLastCommit(void)
{
  // Check if failed at last step
  if (Cfailed)
    return 0;
  else
    return theMaterial->revertToLastCommit();
}

int 
DuctileFracture::revertToStart(void)
{

  Cfailed = false;
  DI		= 0; // Damage index
  DI_VGM	= 0; // Void growth damage component
  DI_MVC	= 0; // Multi-void coalescence damage compoent
  ep_prev	= 0; // Previous plastic strain
  ep_curr	= 0; // Current plastic strain
  dep		= 0; // Incremental plastic strain
  cep_comp	= 0; // Cumulative compressive plastic strain
  es_local = 0; // Local strain
  T = 0; // Triaxiality
  es_max = 0; // The maximum steel strain
  es_min = 0; // The minimum steel strain
  e_memo = 0; // The strain memory factor
  
  Dmax			= 1.0;
  c_mono		= 0;  
  c_cycl		= 0;
  c_symm		= 0;
  E_s			= 0;
  minStrain		= NEG_INF_STRAIN;
  maxStrain		= POS_INF_STRAIN;
  esu = 0;
  k1 = 0;
  k2 = 0;
  db = 0;
  b1 = 0;
  b2 = 0;
  c_dete = 0;

  return theMaterial->revertToStart();
}

UniaxialMaterial *
DuctileFracture::getCopy(void)
{
  DuctileFracture *theCopy = 
    new DuctileFracture(this->getTag(), *theMaterial, Dmax, c_mono, c_cycl, c_symm, E_s, minStrain, maxStrain, esu, k1, k2, db, b1, b2, c_dete);

  theCopy->Cfailed = Cfailed;
  theCopy->trialStrain = trialStrain;

  return theCopy;
}

int 
DuctileFracture::sendSelf(int cTag, Channel &theChannel)
{
    int dbTag = this->getDbTag();

  static ID dataID(3);
  dataID(0) = this->getTag();
  dataID(1) = theMaterial->getClassTag();
  int matDbTag = theMaterial->getDbTag();
  if ( matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    theMaterial->setDbTag(matDbTag);
  }
  dataID(2) = matDbTag;
  if (theChannel.sendID(dbTag, cTag, dataID) < 0) {
    opserr << "DuctileFracture::sendSelf() - failed to send the ID\n";
    return -1;
  }

  static Vector dataVec(25);
  dataVec(0)  = DI;
  dataVec(1)  = DI_VGM;
  dataVec(2)  = DI_MVC;
  dataVec(3)  = ep_prev;
  dataVec(4)  = ep_curr;
  dataVec(5)  = dep;
  dataVec(6)  = cep_comp;
  dataVec(7)  = Dmax;
  dataVec(8)  = c_mono;
  dataVec(9)  = c_cycl;
  dataVec(10) = c_symm;
  dataVec(11) = E_s;
  dataVec(12) = minStrain;
  dataVec(13) = maxStrain;
  dataVec(14) = esu;
  dataVec(15) = k1;
  dataVec(16) = k2;
  dataVec(17) = db;
  dataVec(18) = b1;
  dataVec(19) = b2;
  dataVec(20) = es_max;
  dataVec(21) = es_min;
  dataVec(22) = e_memo;
  dataVec(23) = c_dete;

  if (Cfailed == true)
    dataVec(24) = 1.0;
  else
    dataVec(24) = 0.0;

  if (theChannel.sendVector(dbTag, cTag, dataVec) < 0) {
    opserr << "DuctileFracture::sendSelf() - failed to send the Vector\n";
    return -2;
  }

  if (theMaterial->sendSelf(cTag, theChannel) < 0) {
    opserr << "DuctileFracture::sendSelf() - failed to send the Material\n";
    return -3;
  }

  return 0;
}

int 
DuctileFracture::recvSelf(int cTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static ID dataID(3);
  if (theChannel.recvID(dbTag, cTag, dataID) < 0) {
    opserr << "DuctileFracture::recvSelf() - failed to get the ID\n";
    return -1;
  }
  this->setTag(int(dataID(0)));

  if (theMaterial == 0) {
    int matClassTag = int(dataID(1));
    theMaterial = theBroker.getNewUniaxialMaterial(matClassTag);
    if (theMaterial == 0) {
      opserr << "DuctileFracture::recvSelf() - failed to create Material with classTag " 
	   << dataID(0) << endln;
      return -2;
    }
  }
  theMaterial->setDbTag(dataID(2));

  static Vector dataVec(25);
  if (theChannel.recvVector(dbTag, cTag, dataVec) < 0) {
    opserr << "DuctileFracture::recvSelf() - failed to get the Vector\n";
    return -3;
  }

  DI		= dataVec(0);
  DI_VGM	= dataVec(1);
  DI_MVC	= dataVec(2);
  ep_prev	= dataVec(3);
  ep_curr	= dataVec(4);
  dep		= dataVec(5);
  cep_comp	= dataVec(6);
  Dmax		= dataVec(7);
  c_mono	= dataVec(8);
  c_cycl	= dataVec(9);
  c_symm	= dataVec(10);
  E_s		= dataVec(11);
  minStrain	= dataVec(12);
  maxStrain	= dataVec(13);
  esu		= dataVec(14);
  k1		= dataVec(15);
  k2		= dataVec(16);
  db		= dataVec(17);
  b1		= dataVec(18);
  b2		= dataVec(19);
  es_max	= dataVec(20);
  es_min	= dataVec(21);
  e_memo	= dataVec(22);
  c_dete	= dataVec(23);

  if (dataVec(24) == 1.0)
    Cfailed = true;
  else
    Cfailed = false;

  if (theMaterial->recvSelf(cTag, theChannel, theBroker) < 0) {
    opserr << "DuctileFracture::recvSelf() - failed to get the Material\n";
    return -4;
  }
  return 0;
}

void 
DuctileFracture::Print(OPS_Stream &s, int flag)
{
	if (flag == 100) {
		s << DI << endln;
	}
	
	if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
		s << "DuctileFracture tag: " << this->getTag() << endln;
		s << "\tMaterial: " << theMaterial->getTag() << endln;
		s << "\tDI: " << DI << " Dmax: " << Dmax << endln;
		s << "\tc_mono: " << c_mono << " c_cycl: " << c_cycl << " c_symm: " << c_symm << endln;
	}
		
	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": \"" << this->getTag() << "\", ";
		s << "\"type\": \"DuctileFracture\", ";
		s << "\"material\": \"" << theMaterial->getTag() << "\", ";
		s << "\"tDI\": " << DI << ", ";
		s << "\"Dmax\": " << Dmax << ", ";
		s << "\"tc_mono\": " << c_mono << ", ";
		s << "\"c_cycl\": " << c_cycl << ", ";
		s << "\"c_symm\": " << c_symm << ", ";
	}
}

Response* 
DuctileFracture::setResponse(const char **argv, int argc, OPS_Stream &theOutput)
{
  if (argc == 0) 
    return 0;

  Response *theResponse = 0;

  theOutput.tag("UniaxialMaterialOutput");
  theOutput.attr("matType", this->getClassType());
  theOutput.attr("matTag", this->getTag());


  // stress
  if (strcmp(argv[0],"stress") == 0) {
    theOutput.tag("ResponseType", "sigma11");
    theResponse =  new MaterialResponse(this, 1, this->getStress());
  }  
  // tangent
  else if (strcmp(argv[0],"tangent") == 0) {
    theOutput.tag("ResponseType", "C11");
    theResponse =  new MaterialResponse(this, 2, this->getTangent());
  }

  // strain
  else if (strcmp(argv[0],"strain") == 0) {
    theOutput.tag("ResponseType", "eps11");
    theResponse =  new MaterialResponse(this, 3, this->getStrain());
  }

  // strain
  else if ((strcmp(argv[0],"stressStrain") == 0) || 
	   (strcmp(argv[0],"stressANDstrain") == 0)) {
    theOutput.tag("ResponseType", "sig11");
    theOutput.tag("ResponseType", "eps11");
    theResponse =  new MaterialResponse(this, 4, Vector(2));
  }

 // damage 
  else if (strcmp(argv[0],"damage") == 0) {
    theResponse =  new MaterialResponse(this, 5, DI);
    theOutput.tag("ResponseType", "DI");
  }   

  else if (strcmp(argv[0],"failure") == 0) {
    int res;
    theResponse =  new MaterialResponse(this, 6, res);
    theOutput.tag("ResponseType", "Failure");
  }

  else if (strcmp(argv[0], "vgm") == 0) {
	  int res;
	  theResponse = new MaterialResponse(this, 7, DI_VGM);
	  theOutput.tag("ResponseType", "DI_VGM");
  }

  else if (strcmp(argv[0], "mvc") == 0) {
	  int res;
	  theResponse = new MaterialResponse(this, 8, DI_MVC);
	  theOutput.tag("ResponseType", "DI_MVC");
  }
  // end add


  theOutput.endTag();
  return theResponse;
}

int 
DuctileFracture::getResponse(int responseID, Information &matInfo)
{
  static Vector stressStrain(2);
  static Vector cyclesAndRange(6);

  switch (responseID) {
  case 1:
    matInfo.setDouble(this->getStress());
    return 0;
    
  case 2:
    matInfo.setDouble(this->getTangent());
    return 0;      
    
  case 3:
    matInfo.setDouble(this->getStrain());
    return 0;      
    
  case 4:
    stressStrain(0) = this->getStress();
    stressStrain(1) = this->getStrain();
    matInfo.setVector(stressStrain);
    return 0;
    
  case 5:
    matInfo.setDouble(DI);
    return 0;      

  case 6:
    if (Cfailed == true) 
      matInfo.setInt(1);
    else
      matInfo.setInt(0);
    return 0;   

  case 7:
	  matInfo.setDouble(DI_VGM);
	  return 0;

  case 8:
	  matInfo.setDouble(DI_MVC);
	  return 0;
    
  default:      
    return -1;

  }
}

bool 
DuctileFracture::hasFailed(void) {
  return Cfailed;
}


