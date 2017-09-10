#include<iostream>
#include<vector>
#include<algorithm>
#include<fstream>
#include<map>
#include<cstdlib>

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
};

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

// structure replica
vector<vector<double> > repeat_lattice(vector<vector<double> >,int,int,int);
vector<vector<double> > repeat_coordinate(vector<vector<double> >, int, int, int, int);
vector<vector<double> > repeat_coordinate(vector<vector<double> >,vector<vector<double> >, int, int, int, int);
vector<vector<double> > repeat_xyz(vector<vector<double> >, vector<vector<double> >, int, int, int);
//

// distance
double comp_dist(vector<double>, vector<double>);
//





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
		coord[i][j] = coord_[i][j];
	}
}

/*
cout << lattice[0][0] << "\t" << lattice[0][1]  << "\t" << lattice[0][2] <<endl;
cout << lattice[1][0] << "\t" << lattice[1][1]  << "\t" << lattice[1][2] <<endl;
cout << lattice[2][0] << "\t" << lattice[2][1] << "\t" << lattice[2][2] <<endl;
*/

/* test function lattice_repeated 
vector<vector<double> > lattice_repeated(3,vector<double>(3));
lattice_repeated = repeat_lattice(lattice,3,3,3);
cout << lattice_repeated[0][0] << "\t" << lattice_repeated[0][1]  << "\t" << lattice_repeated[0][2] <<endl;
cout << lattice_repeated[1][0] << "\t" << lattice_repeated[1][1]  << "\t" << lattice_repeated[1][2] <<endl;
cout << lattice_repeated[2][0] << "\t" << lattice_repeated[2][1] << "\t" << lattice_repeated[2][2] <<endl;
*/
// cout << coord[0][0] << "\t" << coord[0][1]  << "\t" << coord[0][2] <<endl;
// cout << coord[1][0] << "\t" << coord[1][1]  << "\t" << coord[1][2] <<endl;
//vector<vector<double> > coords_repeated;
//vector<vector<double> > lattice_repeated(3,vector<double>(3));
vector<vector<double> > xyz_repeated;
//lattice_repeated = repeat_lattice(lattice,3,3,3);
// xyz_repeated= repeat_coordinate(lattice,coord,3,3,3,0);
xyz_repeated= repeat_xyz(lattice,coord,3,3,3);
//xyz_repeated = frac2xyz(lattice_repeated,coords_repeated);
for (int i = 0; i<27; i++){
cout << xyz_repeated[i][0] << "\t" << xyz_repeated[i][1]  << "\t" << xyz_repeated[i][2] <<endl;
//cout << coords_repeated[i][0] << "\t" << coords_repeated[i][1]  << "\t" << coords_repeated[i][2] <<endl;
//cout << coords_repeated[i][0] << "\t" << coords_repeated[i][1] << "\t" << coords_repeated[i][2] <<endl;
}



// end of test block


// vector<vector<double> > lattice(lattice_);
// vector<vector<double> > coord(coord_);
 
//vasp_format(lattice,coord,"BN","POSCAR","direct");
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
vector<vector<double> > ret(fraction.size(),vector<double>(3));
for (int i=0; i<fraction.size();i++){
	ret[i][0]= fraction[i][0]*lattice[0][0] + fraction[i][1] * lattice[1][0] + fraction[i][2] * lattice[2][0];
	ret[i][1]= fraction[i][0]*lattice[0][1] + fraction[i][1] * lattice[1][1] + fraction[i][2] * lattice[2][1];
	ret[i][2]= fraction[i][0]*lattice[0][2] + fraction[i][1] * lattice[1][2] + fraction[i][2] * lattice[2][2];
}
return ret;
}

// unfinished vasp file reading
/*
map<string,vector<vector<double> > > vasp_read(string const file_name){
	ifstream vasp_pos;
	vasp_pos.open(file_name);
	string line;
	vasp_pos
	...
} */

// unfinished vasp file reading

// function to repeat lattice

vector<vector<double> > repeat_lattice(vector<vector<double> > lattice,int nx,int ny,int nz){
	/*
	ret: lattice vectors after repeat;
	lattice: original lattice;
	nx, ny, nz: replicas along three different lattice direction
	*/
	vector<vector<double> > ret(3,vector<double>(3));
	vector<int> repeat(3);
	repeat[0] = nx;
	repeat[1] = ny;
	repeat[2] = nz;
	for (int i=0; i<3; i++){
		for (int j=0; j<3; j++){
			ret[i][j] = lattice[i][j]*repeat[i];
		}
	}
	return ret;
}

// end of function to repeat lattice

// function to repeat coordinates

vector<vector<double> > repeat_coordinate(vector<vector<double> > coords, int nx, int ny, int nz, int type){
	// for output: type = 1 for fractional; type = 0 for xyz;
	// create a replica matrix
	// (1,0,0) (0,1,0), (0,0,1)
	// (1,1,0) (1,0,1), (0,0,1)
	// ...
	int tot_atm = coords.size() * nx * ny * nz;
	vector<vector<double> > ret(tot_atm,vector<double>(3));
	if (type!=1){
		cout << "Please input fractional coordinate!" << endl;
		return ret;
	}
	//cout << tot_atm <<endl;
	for (int i=0; i < nx; i++)
	{
		for (int j=0; j < ny; j++)
		{
			for (int k=0; k<nz; k++)
			{
				for (int atm_i=0; atm_i < coords.size(); atm_i++)
				{
					ret[(i*ny*nz+j*nz+k)*coords.size()+atm_i][0] = (coords[atm_i][0]+i)/nx;
					ret[(i*ny*nz+j*nz+k)*coords.size()+atm_i][1] = (coords[atm_i][1]+j)/ny;
					ret[(i*ny*nz+j*nz+k)*coords.size()+atm_i][2] = (coords[atm_i][2]+k)/nz;
				}
			
			}
		
		}
	
	}
	return ret;

}


// this function can output both fractional and xyz file
// type = 0: xyz coordinate
// type = 1: fractional coordinate
vector<vector<double> > repeat_coordinate(vector<vector<double> > lat,vector<vector<double> > coords, int nx, int ny, int nz,int type){
	// for output: type = 1 for fractional; type = 0 for xyz;
	// create a replica matrix
	// (1,0,0) (0,1,0), (0,0,1)
	// (1,1,0) (1,0,1), (0,0,1)
	// ...
	int tot_atm = coords.size() * nx * ny * nz;
	//cout << tot_atm <<endl;
	vector<vector<double> > ret(tot_atm,vector<double>(3));
	for (int i=0; i < nx; i++)
	{
		for (int j=0; j < ny; j++)
		{
			for (int k=0; k<nz; k++)
			{
				for (int atm_i=0; atm_i < coords.size(); atm_i++)
				{
					ret[(i*ny*nz+j*nz+k)*coords.size()+atm_i][0] = (coords[atm_i][0]+i)/nx;
					ret[(i*ny*nz+j*nz+k)*coords.size()+atm_i][1] = (coords[atm_i][1]+j)/ny;
					ret[(i*ny*nz+j*nz+k)*coords.size()+atm_i][2] = (coords[atm_i][2]+k)/nz;
				}
			
			}
		
		}
	
	}
	if (type==1){
		cout << "Output fractional coordinate!" << endl;
		return ret;
	}
	else {
		cout << "Output xyz coordinate!" << endl;
		vector<vector<double> > ret_xyz;
		vector<vector<double> > ret_lat;
		ret_lat = repeat_lattice(lat,nx,ny,nz);
		ret_xyz = frac2xyz(ret_lat,ret);
		return ret_xyz;
	}
}

// end of function to repeat coordinates

// function to repeat fraction and output xyz file
/*
vector<vector<double> > repeat_coordinate_xyz(vector<vector<double> > lattice,vector<vector<double> > coords, int nx, int ny, int nz){
	int tot_atm = coords.size() * nx * ny * nz;
	int type = 1;
	//cout << tot_atm <<endl;
	vector<vector<double> > lat(3,vector<double>(3));
	vector<vector<double> > coord(tot_atm,vector<double>(3));
	vector<vector<double> > ret_xyz(tot_atm,vector<double>(3));
	lat=repeat_lattice(lattice,nx,ny,nz);
	coord=repeat_coordinate(coords,nx,ny,nz,type);
	ret_xyz = frac2xyz(lat,coord);
	return ret_xyz;
}
*/

// function to repeat xyz file
vector<vector<double> > repeat_xyz(vector<vector<double> > lat, vector<vector<double> > coords, int nx, int ny, int nz)
{
	int tot_atm = coords.size() * nx * ny * nz;
	//cout << tot_atm <<endl;
	vector<vector<double> > ret(tot_atm,vector<double>(3));
	for (int i=0; i < nx; i++)
	{
		for (int j=0; j < ny; j++)
		{
			for (int k=0; k<nz; k++)
			{
				for (int atm_i=0; atm_i < coords.size(); atm_i++)
				{
					ret[(i*ny*nz+j*nz+k)*coords.size()+atm_i][0] = coords[atm_i][0]+i*lat[0][0]+j*lat[1][0]+k*lat[2][0];
					ret[(i*ny*nz+j*nz+k)*coords.size()+atm_i][1] = coords[atm_i][1]+i*lat[0][1]+j*lat[1][1]+k*lat[2][1];
					ret[(i*ny*nz+j*nz+k)*coords.size()+atm_i][2] = coords[atm_i][2]+i*lat[0][2]+j*lat[1][2]+k*lat[2][2];
				}
			
			}
		
		}
	
	}
	return ret;
}

// end function

/* function to calculate distance between atoms

 input: lat: lattice vector of the system
 		type: whether the system is periodic or not
 		dims: dims = 0 indicate a 0D system, no periodic condition
 			  dims = 1 1D system
 			  dims = 2 2D system
 			  dims = 3 3D system
 return: a n_atom by n_atom matrix which store the distance between atomic pairs
*/

double comp_dist(vector<double> pos1, vector<double> pos2){
	return sqrt((pos1[0]-pos2[0])**2+(pos1[1]-pos2[1])**2+(pos1[2]-pos2[2])**2);
}

vector<vector<double> > dist_map (vector<vector<double> > lat, vector<vector<double> > coords, int dims){
	
	vector<vector<double> > ret_dist(coords.size(),vector<double>(coords.size(),0.0));
	
	if (dims == 0)
	{
		for (int i=0; i<coords.size()-1; i++)
		{
			for (int j=i+1; j< coords.size(); j++)
			{
				ret_dist[i][j]=comp_dist(coords[i],coords[j]);
				ret_dist[j][i]= ret_dist[i][j];
			}
		}
		return ret_dist;
	}
	else if (dims == 1)
	{
		// assume along a1 direction
		// generate a |_|c|_| box
		vector<vector<double> > rep_m<coords.size(),vector<double>(3)>; // -a1
		vector<vector<double> > rep_p<coords.size(),vector<double>(3)>; // +a1
		for (int i=0; i< coords.size(); i++)
		{
			rep_m[i][0]= coords[i][0]-lat[0][0];
			rep_m[i][1]= coords[i][1]-lat[0][1];
			rep_m[i][2]= coords[i][2]-lat[0][2];
			rep_p[i][0]= coords[i][0]+lat[0][0];
			rep_p[i][1]= coords[i][1]+lat[0][1];
			rep_p[i][2]= coords[i][2]+lat[0][2];
			
		}
		for (int i=0; i< coords.size()-1; i++)
		{
			for (int j=i+1; j<coords.size(); j++)
			{
				double temp = comp_dist(coords[i],coords[j]);
				temp = temp<comp_dist(coords[i],rep_m[j])?temp:comp_dist(coords[i],rep_m[j]);
				temp = temp<comp_dist(coords[i],rep_p[j])?temp:comp_dist(coords[i],rep_p[j]);
				ret_dist[i][j]=temp;
				ret_dist[j][i]=temp;
			}
		}
		return ret_dist;
	}
	else if (dims == 2)
	{
		// assume repeate along a1 and a2 direction
		// we set up a 9 to represent the neighboring box
		// matrix[1][1] is our center box or original structure
		vector<vector<vector<double> > > square(9,vector<vector<double> >(coords.size(),vector<double>(3)));
		int a1_op[9] = {0,0,0,1,1,1,-1,-1,-1};
		int a2_op[9] = {0,1,-1,0,1,-1,0,1,-1};
		for (int i=0; i<9; i++){
			for (int j=0; j<coords.size(); j++)
			{
				for (int k=0; k<3; k++)
				{
					square[i][j][k] = coords[j][k] + a1_op[i]*lat[0][k] +a2_op[i]*lat[1][k];
				}
			}
		}
		for (int i = 0; i<coords.size()-1; i++)
		{
			for (int j=i; j< coords.size(); j++)
			{
				double temp_min = comp_dist(coords[i],coords[j]);
				for (int k=1; k<9; k++)
				{
					temp_min = min(temp_min,comp_dist(coords[i],square[k][j]));
				}
				ret_dist[i][j]=temp_min;
				ret_dist[j][i]=ret_dist[i][j];
			}
		}
		return ret_dist;
	}
	else if (dims == 3)
	{
		// assume repeate along a1, a2, and a3 direction
		// we set up a 27 to represent the neighboring box
		// matrix[1][1] is our center box or original structure
		vector<vector<vector<double> > > cubic(27,vector<vector<double> >(coords.size(),vector<double>(3)));
		int a1_op[27] = {0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
		int a2_op[27] = {0,0,0,1,1,1,-1,-1,-1,0,0,0,1,1,1,-1,-1,-1,0,0,0,1,1,1,-1,-1,-1};
		int a3_op[27] = {0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1};
		for (int i=0; i<27; i++){
			for (int j=0; j<coords.size(); j++)
			{
				for (int k=0; k<3; k++)
				{
					cubic[i][j][k] = coords[j][k] + a1_op[i]*lat[0][k] +a2_op[i]*lat[1][k] + a3_op[i]*lat[2][k];
				}
			}
		}
		for (int i = 0; i<coords.size()-1; i++)
		{
			for (int j=i; j< coords.size(); j++)
			{
				double temp_min = comp_dist(coords[i],coords[j]);
				for (int k=1; k<27; k++)
				{
					temp_min = min(temp_min,comp_dist(coords[i],cubic[k][j]));
				}
				ret_dist[i][j]=temp_min;
				ret_dist[j][i]=ret_dist[i][j];
			}
		}
		return ret_dist;
	}
	else
	{
		cout << "Dimension > 3" << endl;
		return ret_dist;
	}
	
}