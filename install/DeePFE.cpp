#include "Bias.h"
#include "ActionRegister.h"

#include "tensorflow/core/public/session.h"
#include "tensorflow/core/platform/env.h"
#include "tensorflow/core/framework/op.h"
#include "tensorflow/core/framework/op_kernel.h"
#include "tensorflow/core/framework/shape_inference.h"

using namespace std;
using namespace tensorflow;

typedef float  VALUETYPE;

namespace PLMD {
  namespace bias {

    const double global_deepfe_f_cvt = 96.485;
    const double global_deepfe_e_cvt = 96.485;

    class DeePFE : public Bias {
  public:
      explicit DeePFE(const ActionOptions&);
      virtual ~DeePFE();
      void calculate();
      static void registerKeywords(Keywords& keys);
  private:
      Session** session;
      GraphDef* graph_def;
      int numb_model;
      vector<string> graph_names;
      void checkStatus(const tensorflow::Status& status) ;
      Value* valueForce2;
      Value* valueFstd;
      void compute_std (vector<double > & avg_f,
			vector<double > & std_f,
			const vector<vector<double > > ff) const;
      double avg_std (const vector<double > & std_f) const;
      double trust_l_1, trust_l_2;
      double deepfe_e_cvt, deepfe_f_cvt;
      double apply_trust (const double & xx) const;
      bool use_scaled_f;
    };

    PLUMED_REGISTER_ACTION(DeePFE,"DEEPFE")

    void DeePFE::checkStatus(const tensorflow::Status& status) {
      if (!status.ok()) {
	std::cout << status.ToString() << std::endl;
	exit(1);
      }
    }

    void DeePFE::registerKeywords(Keywords& keys) {
      Bias::registerKeywords(keys);
      keys.use("ARG");
      keys.addFlag("SCALED_FORCE", false, "whether to use scaled force to calculate std");
      keys.add("compulsory","MODEL","graph.pb","the model file of DeePFE");
      keys.add("compulsory","TRUST_LVL_1","the trust level 1 of the force");
      keys.add("compulsory","TRUST_LVL_2","the trust level 2 of the force");
      keys.add("optional","UNIT_CVT","The unit conversion constant of energy from graph to kJ/mol (default is 96.485: from eV to kJ/mol)");
      keys.addOutputComponent("force2","default","the instantaneous value of the squared force due to this bias potential");
      keys.addOutputComponent("fstd","default","the standard deviation of the biasing force from the ensemble of models");
    }

    DeePFE::DeePFE (const ActionOptions&ao)
	: PLUMED_BIAS_INIT(ao), 
	  deepfe_e_cvt(global_deepfe_e_cvt), 
	  deepfe_f_cvt(global_deepfe_f_cvt)
    {
      parseFlag("SCALED_FORCE", use_scaled_f);
      parseVector ("MODEL", graph_names);
      parse ("TRUST_LVL_1", trust_l_1);
      parse ("TRUST_LVL_2", trust_l_2);
      assert (trust_l_1 <= trust_l_2);
      parse("UNIT_CVT", deepfe_f_cvt);
      deepfe_e_cvt = deepfe_f_cvt;
      checkRead();
      numb_model = graph_names.size();

      log.printf("  graphname");
      for (int ii = 0; ii < numb_model; ++ii){
	      log.printf(" %s", graph_names[ii].c_str());
      }
      log.printf("\n");
      log.printf("  trust levels: %f %f\n", trust_l_1, trust_l_2);
      log.printf("  unit_cvt: %f\n", deepfe_e_cvt);

      addComponent("force2");
      componentIsNotPeriodic("force2");
      valueForce2=getPntrToComponent("force2");

      addComponent("fstd");
      componentIsNotPeriodic("fstd");
      valueFstd=getPntrToComponent("fstd");

      session = new Session*[numb_model];
      graph_def = new GraphDef[numb_model];
      for (int ii = 0; ii < numb_model; ++ii){
	checkStatus (NewSession(SessionOptions(), &session[ii]));
	checkStatus (ReadBinaryProto(Env::Default(), graph_names[ii], &graph_def[ii]));
	checkStatus (session[ii]->Create(graph_def[ii]));
      }
    }

    DeePFE::~DeePFE ()
    {
      delete [] session;
      delete [] graph_def;
    }

    void DeePFE::calculate() {
      double ene=0.0;
      double totf2=0.0;
      vector<double > xx(getNumberOfArguments());
      vector<vector<double > > dforce (numb_model);
      vector<double> dener (numb_model, 0);
      vector<vector<double > > scaledforce (numb_model); 

      for(unsigned ii=0; ii<getNumberOfArguments(); ++ii) {
	xx[ii] = getArgument(ii);
      }

      {
	int nframes = 1;
	int isize = xx.size();

	TensorShape coord_shape ;
	coord_shape.AddDim (nframes);
	coord_shape.AddDim (isize * 2);  

	Tensor coord_tensor	(DT_FLOAT, coord_shape);

	auto coord = coord_tensor.matrix<VALUETYPE> ();

	for (int ii = 0; ii < nframes; ++ii){
	  for (int jj = 0; jj < isize; ++jj){
	    coord(ii, jj) = xx[jj];
	  }
	  for (int jj = isize; jj < isize * 2; ++jj){
	    coord(ii, jj) = 0.;
	  }
	}

  Tensor drop_out_rate_tensor (DT_FLOAT, TensorShape({1}));
  auto drop_out_rate = drop_out_rate_tensor.scalar<VALUETYPE> ();
  drop_out_rate(0) = 0.;

	std::vector<std::pair<string, Tensor>> input_tensors = {
	  {"inputs",	coord_tensor}, 
    {"drop_out_rate", drop_out_rate_tensor}
	};
	std::vector<Tensor> output_tensors;

	for (int kk = 0; kk < numb_model; ++kk) {
	  checkStatus (session[kk]->Run(input_tensors, 
					{"o_energy", "o_forces"}, 
					{}, 
					&output_tensors));
  
	  Tensor output_e = output_tensors[0];
	  Tensor output_f = output_tensors[1];
	  auto oe = output_e.flat <VALUETYPE> ();
	  auto of = output_f.flat <VALUETYPE> ();
	  dener[kk] = oe(0) * deepfe_e_cvt;
	  dforce[kk].resize(isize);
	  for (int ii = 0; ii < isize; ++ii){
	    dforce[kk][ii] = of(ii) * deepfe_f_cvt;
	  }
          if (use_scaled_f) {
            std::vector<Tensor> output_scaled_f;
            checkStatus (session[kk]->Run(input_tensors,
                                          {"o_scaled_forces"},
                                          {},
                                          &output_scaled_f));
            auto osf = output_scaled_f[0].flat <VALUETYPE> ();
            scaledforce[kk].resize(isize);
            for (int ii = 0; ii < isize; ++ii){
              scaledforce[kk][ii] = osf(ii);
            }
          }
	  // for (int ii = 0; ii < isize; ++ii){
	  //   cout << xx[ii] << " " ;
	  // }
	  // cout << "    " ;
	  // cout << dener[kk] ;
	  // cout << "    " ;
	  // for (int ii = 0; ii < isize; ++ii){
	  //   cout << dforce[kk][ii] << " " ;
	  // }
	  // cout << endl;
	}
      }
      // average energy
      double avg_e = 0;
      for (unsigned ii = 0; ii < dener.size(); ++ii){
	avg_e += dener[ii];
      }
      avg_e /= double(dener.size());
      // average and std of force
      vector<double > avg_f, std_f;
      compute_std (avg_f, std_f, dforce);
      double avg_std_f = avg_std (std_f);
      if (use_scaled_f) {
        vector<double > avg_sf, std_sf;
        compute_std (avg_sf, std_sf, scaledforce);
        avg_std_f = avg_std (std_sf);
      }
      // compute trust
      double trust = apply_trust (avg_std_f);
      // for (unsigned ii = 0; ii < xx.size(); ++ii){
      // 	cout << xx[ii] << " " ;
      // }
      // cout << "  \t" << avg_std_f 
      // 	   << "  \t" << trust
      // 	   << endl;

      ene+=trust * avg_e;
      // cout << "real e "  << trust * avg_e << endl;
      for(unsigned ii=0; ii<getNumberOfArguments(); ++ii) {
	const double f=trust * avg_f[ii];
	// cout << "trust " << trust << " orig f " << avg_f[ii] << " real f "  << f << endl;
	setOutputForce(ii,-f);
	totf2+=f*f;
      }
      setBias(-ene);
      valueForce2->set(totf2);
      valueFstd->set(avg_std_f);
    }

    void DeePFE::compute_std (
	vector<double > & avg_f,
	vector<double > & std_f,
	const vector<vector<double > > ff) const
    {
      int numb_comp = ff[0].size();
      avg_f.clear();
      avg_f.resize (numb_comp, 0);
      assert (ff.size() == numb_model);
      for (int ii = 0; ii < numb_model; ++ii){
	for (unsigned jj = 0; jj < ff[ii].size(); ++jj){
	  avg_f[jj] += ff[ii][jj];
	}
      }
      for (unsigned ii = 0; ii < avg_f.size(); ++ii){
	avg_f[ii] /= double (numb_model);
	// cout << avg_f[ii] << " " ;
      }
      // cout << endl;
      std_f.clear();
      std_f.resize (numb_comp, 0);
      for (int ii = 0; ii < numb_model; ++ii){
	for (int jj = 0; jj < numb_comp; ++jj){
	  std_f[jj] += (ff[ii][jj] - avg_f[jj]) * (ff[ii][jj] - avg_f[jj]);
	}
      }      
      // double avg_std = 0;
      for (unsigned ii = 0; ii < std_f.size(); ++ii){
	std_f[ii] /= double (numb_model);
	std_f[ii] = sqrt(std_f[ii]);
	// avg_std += std_f[ii];
	// cout << std_f[ii] << " " ;
      }
      // cout << endl;
      // avg_std /= double(std_f.size());
      // cout << avg_std << endl;
    }

    double DeePFE::avg_std (const vector<double > & std_f) const
    {
      double avg = 0;
      for (unsigned ii = 0; ii < std_f.size(); ++ii){
	avg += std_f[ii] * std_f[ii];
      }
      return sqrt (avg / double(std_f.size()));
    }

    double DeePFE::apply_trust (const double & xx) const
    {
      const double &rmin = trust_l_1;
      const double &rmax = trust_l_2;
      double vv = 0;
      if (xx < rmin) {
	vv = 1;
      }
      else if (xx < rmax){
	double value = (xx - rmin) / (rmax - rmin) * M_PI;
	vv = 0.5 * (cos(value) + 1);
      }
      else {
	vv = 0;
      }
      return vv;
    }
  }
}

