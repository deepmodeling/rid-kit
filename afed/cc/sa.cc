#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <queue>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;
#define MaxLineLength 65536

namespace StringOperation{
  void split (const std::string & in,
	      std::vector<std::string > & out);
}


void StringOperation::
split (const std::string & in,
       std::vector<std::string > & out)
{
  std::istringstream iss(in);
  out.clear();
  
  do {
    std::string sub;
    iss >> sub;
    out.push_back (sub);
  } while (iss);

  out.pop_back();
}

void diff_pbc (vector<double > & diff, 
	       const vector<double > & aa,
	       const vector<double > & bb) 
{
  unsigned nvalue = aa.size();
  diff.resize (nvalue);
  for (unsigned ii = 0; ii < nvalue; ++ii){
    diff[ii] = aa[ii] - bb[ii];
    if (diff[ii] > 180) {
      diff[ii] -= 360;
    }
    else if (diff[ii] < -180) {
      diff[ii] += 360;
    }
  }
}

double norm (const vector<double > & vv)
{
  double ret = 0;
  for (unsigned ii = 0; ii < vv.size(); ++ii){
    ret += vv[ii] * vv[ii];
  }
  return sqrt(ret);
}

void load_table (vector<vector<double > > & odata,
		 const string & ifile) 
{
  ifstream data ;
  data.open (ifile.c_str());
  if (!data){
    cerr << "cannot open file \"" << ifile << "\"" << endl;
    throw;
  }

  odata.clear();
  char valueline [MaxLineLength];
  while (data.getline(valueline, MaxLineLength)){
    if (valueline[0] == '#' || valueline[0] == '@'){
      continue;
    }
    vector<string > words;
    StringOperation::split (string(valueline), words);
    vector<double > vline (words.size());
    for (unsigned ii = 0; ii < words.size(); ++ii){
      vline[ii] = (atof(words[ii].c_str()));
    }
    odata.push_back(vline);
  }
}


int main(int argc, char * argv[])
{
  if (argc != 2) {
    cerr << "usage" << endl;
    cerr << argv[0] << " threshold" << endl;
    return 1;
  }
  double threshold = atof(argv[1]);  

  string ifile = "angle.deg.out";
  string ofile = "sel.out";
  string oifile = ofile + ".idx";
  
  cout << "input:     " << ifile << endl;
  cout << "output:    " << ofile << " " << oifile << endl;
  cout << "threshold: " << threshold << endl;

  vector<vector<double > > data;
  load_table (data, ifile);
  
  vector<int > sel_idx;
  sel_idx.push_back (0);
  
  for (unsigned ii = 1; ii < data.size(); ++ii) {
    bool do_append = true;
    for (unsigned jj = 0; jj < sel_idx.size(); ++jj){
      vector<double > diff;
      diff_pbc (diff, data[ii], data[sel_idx[jj]]);
      if (norm(diff) < threshold) {
	do_append = false;
	break;
      }
    }
    if (do_append){
      sel_idx.push_back (ii);
      cout << ii << " " << sel_idx.size() << endl;
    }
  }

  ofstream ofp (ofile);
  ofstream ofpi (oifile);
  
  for (unsigned ii = 0; ii < sel_idx.size(); ++ii){
    for (unsigned jj = 0; jj < data[sel_idx[ii]].size(); ++jj){
      ofp << data[sel_idx[ii]][jj] << " " ;
    }
    ofp << endl;
    ofpi << sel_idx[ii] << endl;
  }

  return 0;
}
