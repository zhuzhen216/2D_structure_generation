#include<iostream>
#include<vector>
#include<algorithm>
#include<fstream>
#include<map>

using namespace std;


// *************************************************
//
// Declarations
//
// *************************************************
/*
x, y, z coordinate structure
*/

struct coordination{
double x;
double y;
double z;
}

/*
vasp_read function: read a POSCAR file
*/

// compound_split would split a string to show how many elements in the
// compound, and give the fraction of the compound
map<string,double> compound_split(string const);


// output input file for VASP
void vasp_format(vector<vector<double> >&,vector<vector<double> >&,string,string,string);


// output input file for SIESTA
void siesta_format(vector<vector<int> >&,vector<vector<int> >&,string,string,string);

map<string,vector<vector<double> > > vasp_read(string const);

//
vector<vector<double> > frac2xyz(vector<vector<double> >, vector<vector<double> >);
vector<coordination> frac2xyz(vector<vector<double> >, vector<coordination>);



// *************************************************
//
// Main block
//
// *************************************************



int main(){


// test block
// ofstream structure;
double lattice_[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
double coord_[2][3]={{0,0,0},{0.5,0.5,0}};

vector<vector<double> > lattice(3,vector<double>(3));
vector<vector<double> > coord(2,vector<double>(3));

for (int i=0; i< lattice.size(); i++){
	for (int j=0;j< lattice[0].size(); j++){
		lattice[i][j] = lattice_[i][j];
	}
}

for (int i=0; i< coord.size(); i++){
	for (int j=0;j< coord[0].size(); j++){
		lattice[i][j] = coord_[i][j];
	}
}
// end of test block


// vector<vector<double> > lattice(lattice_);
// vector<vector<double> > coord(coord_);
 
vasp_format(lattice,coord,"BN","POSCAR","direct");
/*
structure.open("test.xyz");
structure<< "test for structure generation" <<endl;
structure.close();
*/
return 0;
}


/*

function to output vasp_format 

*/


void vasp_format(vector<vector<double> >& lattice,vector<vector<double> >& coord,string compound, string name="POSCAR",string type="direct"){
/* 
Inputs:
lattice: the lattice vector of the structure
coord: fractional coordinate of the structure by default
compound: the compound name
name: the file name
type: fractional or card

Variables:
vasp_pos: the ofstream file;
natm: number of atoms in the file;
compound_info: take BN for example {"B":0.5;"N":0.5}

*/
// default to pass a fractional coordinate
// compound name follow the case: P, AlB,SiS, and so on
// try P first
ofstream vasp_pos;
// scale : whether apply scaling on the POSCAR
double scale = 1.0;
//nat: number of atoms in the POSCAR
int nat = coord.size();
vasp_pos.open(name);
vasp_pos << compound << endl;
vasp_pos << scale << endl;
for (int i=0; i<lattice.size(); i++){
	for(int j=0; j<lattice[0].size();j++){
		vasp_pos << lattice[i][j] << '\t';
}
	vasp_pos << endl;
}

map<string,double> compound_info;
compound_info = compound_split(compound);
for(map<string,double>::const_iterator iter=compound_info.begin();iter!=compound_info.end();iter++){
	vasp_pos << iter->first << '\t';
}
vasp_pos << endl;

for(map<string,double>::const_iterator iter=compound_info.begin();iter!=compound_info.end();iter++){
	vasp_pos << nat*iter->second << '\t';
}

vasp_pos << endl;
vasp_pos << type <<endl;
for (int i = 0; i < coord.size();i++){
	for (int j=0; j<3;j++){
		vasp_pos << coord[i][j] << '\t';
}
vasp_pos<<endl;
}

vasp_pos.close();
}


/*****

function to split a compound name and extract element information

****/
map<string,double> compound_split(string const compound){
	// input P, output {P:1.0}
	map<string,double> ret;

	// case 1: for elemental materials
	if (compound.length()==1){
		ret[compound]=1.0;
		return ret;
	}

	int i = 0;
	int count = 0; // count how many elements in the compound
	int sum_atom = 0;
	while (i<compound.length())
	{
		if (isupper(compound[i])){
			int start = i;
			int length = 1;
			i++;
			while (i<compound.length()&&islower(compound[i])){
				i++;
				length++;
			}
			string elem = compound.substr(start,length);
			if (i==compound.length() || isupper(compound[i])){
				ret[elem]=1;
				sum_atom+=ret[elem];
			}
			else{
				int digit_start = i;
				int digit_length = 1;
				while (i<compound.length()&&isdigit(compound[i])){
					i++;
					digit_length++;
				}
				ret[elem]=stoi(compound.substr(digit_start,digit_length));
				sum_atom+=ret[elem];
			}
	}
	}
    //*	
	for (map<string,double>::const_iterator iter=ret.begin();iter!=ret.end();iter++){
		ret[iter->first] = ret[iter->first]/sum_atom;
	}
    //*/
	return ret;
}

/*
function: convert fractional coordinate to xyz using coordinate structure
*/

vector<coordination> frac2xyz(vector<vector<double> > lattice, vector<coordination> fraction){
vector<coordination> ret(fraction.size());
for (int i=0; i<fraction.size();i++){
	ret[i].x= fraction[i].x*lattice[0][0] + fraction[i].y * lattice[1][0] + fraction[i].z * lattice[2][0];
	ret[i].y= fraction[i].x*lattice[0][1] + fraction[i].y * lattice[1][1] + fraction[i].z * lattice[2][1];
	ret[i].z= fraction[i].x*lattice[0][2] + fraction[i].y * lattice[1][2] + fraction[i].z * lattice[2][2];
}
return ret;
}

/*
function: convert fractional coordinate to xyz in the form of vectors
*/

vector<vector<double> > frac2xyz(vector<vector<double> > lattice, vector<vector<double> > fraction){
vector<coordination> ret(fraction.size(),vector<double>(3));
for (int i=0; i<fraction.size();i++){
	ret[i][0]= fraction[i][0]*lattice[0][0] + fraction[i][1] * lattice[1][0] + fraction[i][2] * lattice[2][0];
	ret[i][1]= fraction[i][0]*lattice[0][1] + fraction[i][1] * lattice[1][1] + fraction[i][2] * lattice[2][1];
	ret[i][2]= fraction[i][0]*lattice[0][2] + fraction[i][1] * lattice[1][2] + fraction[i][2] * lattice[2][2];
}
return ret;
}


map<string,vector<vector<double> > > vasp_read(string const file_name){
	ifstream vasp_pos;
	vasp_pos.open(file_name);
	string line;
	vasp_pos
}

