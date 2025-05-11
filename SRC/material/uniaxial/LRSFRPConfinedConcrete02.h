#ifndef LRSFRPConfinedConcrete02_h
#define LRSFRPConfinedConcrete02_h

#include <UniaxialMaterial.h>

class LRSFRPConfinedConcrete02 : public UniaxialMaterial
{
public:
	LRSFRPConfinedConcrete02(int tag, double fc0, double Ec, double ec0, double t, double Efrp1, double Efrp2, double eps_frp0, double eps_h_rup, double C, double R, double ft, double Ets, int Unit);
	LRSFRPConfinedConcrete02(int tag, double fc0, double Ec, double ec0, double fcc, double ecu, double ft, double Ets, int Unit);
	LRSFRPConfinedConcrete02(int tag, double fc0, double Ec, double ec0, double ft, double Ets, int Unit);
	LRSFRPConfinedConcrete02();

	~LRSFRPConfinedConcrete02();

	const char* getClassType(void) const { return "LRSFRPConfinedConcrete02"; }

	int setTrialStrain(double strain, double strainRate = 0.0);
	double getStrain(void);
	double getStress(void);
	double getTangent(void);

	double getInitialTangent(void) { return m_Ec; };

	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	UniaxialMaterial* getCopy(void);

	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel,
		FEM_ObjectBroker& theBroker);

	void Print(OPS_Stream& s, int flag = 0);

	//////////////////////////////////////////////////////////////////////////
	void Tens_Envlp(double epsc, double& sigc, double& Ect);
	void Compr_Envlp(double epsc, double& sigc, double& Ect);
	void Compr_UnloadingPath(double epsc, double& sigc, double& Ect);
	void Compr_GetPlasticStrain();
	void Compr_GetStrainRecoveryRatio();
	void Compr_ReloadingPath(double epsc, double& sigc, double& Ect);
	void GetRefPoint();
	void GetDeterioratedStress();
	void GetStressDeteriorationRatio();

	// AddingSensitivity:BEGIN //////////////////////////////////////////
	int    setParameter(const char** argv, int argc, Information& info);
	int    updateParameter(int parameterID, Information& info);
	int    activateParameter(int parameterID);
	double getStressSensitivity(int gradNumber, bool conditional);
	int    commitSensitivity(double strainGradient, int gradNumber, int numGrads);
	// AddingSensitivity:END ///////////////////////////////////////////
protected:

private:
	//////////////////////////////////////////////////////////////////////////
	//inputs
	double m_fc0;
	double m_Ec;
	double m_epsc0;
	double m_t;
	double m_Efrp1;
	double m_Efrp2;
	double m_eps_frp0;
	double m_eps_h_rup;
	double m_C;
	double m_R;
	double m_Ets;
	double m_ft;

	int m_Unit;
	//////////////////////////////////////////////////////////////////////////
	//Assistant variable
	double m_rou;
	double m_rou2;
	double m_epst;
	double m_epst2;
	double m_epscu;
	double m_fl;
	double m_fcc;
	double m_E2;
	double m_E22;
	double m_Unitscale;
	//////////////////////////////////////////////////////////////////////////
	//Current state variable
	//Tension
	double m_epstn;
	double m_epstu;
	double m_Etr1;
	double m_Etr2;

	//Compression
	int m_n;
	int m_ne;
	int m_loadingflag;

	double m_Ere;
	double m_epsunenv;
	double m_sigunenv;
	double m_sigunn1;
	double m_epsretenv;
	double m_epsun;
	double m_sigun;
	double m_epsre;
	double m_sigre;
	double m_epsref;
	double m_sigref;
	double m_Eun0;

	double m_betaun;
	double m_fi;
	double m_fiful;
	double m_signew;

	double m_gamare;
	double m_omg;
	double m_omgful;
	double m_epspl;

	bool m_bSmallStress;
	bool m_bUltimate;

	double m_Tstrain;
	double m_Tstress;

	double m_trialTangent;

	//////////////////////////////////////////////////////////////////////////
	//History variable
	//Tension
	double m_epstnlast;
	double m_epstulast;
	double m_Etr1last;
	double m_Etr2last;

	//Compression
	int m_nlast;
	int m_nelast;
	int m_loadingflaglast;

	double m_Erelast;
	double m_epsunenvlast;
	double m_sigunenvlast;
	double m_sigunn1last;
	double m_epsretenvlast;
	double m_epsunlast;
	double m_sigunlast;
	double m_epsrelast;
	double m_sigrelast;
	double m_epsreflast;
	double m_sigreflast;
	double m_Eun0last;

	double m_betaunlast;
	double m_filast;
	double m_fifullast;
	double m_signewlast;

	double m_gamarelast;
	double m_omglast;
	double m_omgfullast;
	double m_epspllast;

	bool m_bSmallStresslast;
	bool m_bUltimatelast;

	double m_trialStrainlast;
	double m_trialStresslast;

	double m_trialTangentlast;

	// AddingSensitivity:BEGIN //////////////////////////////////////////
	int parameterID;
	Matrix* SHVs;
	// AddingSensitivity:END ///////////////////////////////////////////
};
#endif
