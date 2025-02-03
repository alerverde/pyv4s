#include "Vector.h"
#include "Tetrahedron.h"
#include <iostream>
#include <algorithm>
#include <cstring>

float** cosort(std::vector<float> values, float xs[4], float ys[4], float zs[4]) {
    int ids[4]= {0,1,2,3};
    std::vector<std::pair<float,int>> pairs;
    for(int i= 0; i < 4; i++)
        pairs.emplace_back(values[i],ids[i]);
    std::sort(pairs.begin(),pairs.end(),[](const std::pair<float,int>& a, const std::pair<float,int>& b) {
        return a.first < b.first;
    });

    float** sorted_list= new float*[4];
    for(int i= 0; i < 4; i++) {
        sorted_list[i]= new float[4];
        int id= pairs[i].second;
        sorted_list[i][0]= xs[id];
        sorted_list[i][1]= ys[id];
        sorted_list[i][2]= zs[id];
        sorted_list[i][3]= pairs[i].first;
    }
    return sorted_list;
}

float getPotential(Atom a1, Atom a2, Vector bounds) {
    float R= distanceBetween(a1,a2,bounds);
    float e= sqrt(a1.e * a2.e);
    float s= .5*(a1.s + a2.s);
    const float K= 1389.35458; //kJ/mol to e-2/A

    float output= 0;
    if(e != 0) output= 4*e*(pow(s/R,12)-pow(s/R,6));
    output+= K * a1.q * a2.q / R;
    return output;
}

Vector parseBounds(float* bounds_f) {
    Vector bounds;
    bounds.x= bounds_f[0];
    bounds.y= bounds_f[1];
    bounds.z= bounds_f[2];
    return bounds;
}

bool shouldSkipAtom(int a, int& neighbor, int atoms_per_water_molecule, const Atom* atom_list, const bool neighbor_is_O, Vector bounds, const float R_CUT_OFF) {
    if(neighbor >= a && neighbor < a + atoms_per_water_molecule)
        return true;

    if(distanceBetween(atom_list[neighbor], atom_list[a], bounds) > R_CUT_OFF+1.1) {
        if(neighbor_is_O)
            neighbor+= atoms_per_water_molecule - 1;
        return true;
    }

    return false;
}

int getIndexOfClosestPoint(int neighbor, const Atom* atom_list, std::vector<Vector> sites, Vector bounds, const float R_CUT_OFF) {
    int index_closest= 0;
    float dist_closest= distanceBetween(atom_list[neighbor],sites[index_closest],bounds);

    for(std::size_t idx_site= 1; idx_site < sites.size(); idx_site++) {
        float dist_site= distanceBetween(atom_list[neighbor],sites[idx_site],bounds);
        if(dist_site < dist_closest) {
            dist_closest= dist_site;
            index_closest= idx_site;
        }
    }

    if(dist_closest > R_CUT_OFF)
        index_closest= -1;
    return index_closest;
}

float getPairwaysInteractions(int a, int& neighbor, bool neighbor_is_O, const Atom* atom_list, int atoms_per_water_molecule, Vector bounds) {
    float interactions= 0.;
    for(int idx_atom= a; idx_atom < a+atoms_per_water_molecule; idx_atom++)
        interactions+= getPotential(atom_list[idx_atom],atom_list[neighbor],bounds);
    
    if(!neighbor_is_O)
        return interactions;

    for(int idx_atom_a= a; idx_atom_a < a+atoms_per_water_molecule; idx_atom_a++)
        for(int idx_atom_neighbor= neighbor+1; idx_atom_neighbor < neighbor+atoms_per_water_molecule; idx_atom_neighbor++)
            interactions+= getPotential(atom_list[idx_atom_a],atom_list[idx_atom_neighbor],bounds);
    neighbor+= atoms_per_water_molecule-1;

    return interactions;
}

std::vector<float> getInteractionPerSite(const int a, Atom* atom_list, std::vector<Vector> sites, const int N_ATOMS, const int atoms_per_water_molecule, const char* type_O, Vector bounds, const float R_CUT_OFF) {
    std::vector<float> sum_site(4);

    for(int neighbor= 0; neighbor < N_ATOMS; neighbor++) {
        bool neighbor_is_O= strcmp(atom_list[neighbor].a, type_O) == 0;
        if(shouldSkipAtom(a, neighbor, atoms_per_water_molecule, atom_list, neighbor_is_O, bounds, R_CUT_OFF)) continue;

        int index_closest= getIndexOfClosestPoint(neighbor, atom_list, sites, bounds, R_CUT_OFF);
        if(index_closest == -1) {
            if(neighbor_is_O)
                neighbor+= atoms_per_water_molecule -1;
            continue;
        }

        sum_site[index_closest]+= getPairwaysInteractions(a, neighbor, neighbor_is_O, atom_list, atoms_per_water_molecule, bounds);
    }

    return sum_site;
}


extern "C" float*** tetrahedrons(Atom* atom_list, int N_ATOMS, float* bounds_f, bool* study_molecules, const char* type_O, const int atoms_per_water_molecule, const float R_CUT_OFF) {
    Vector bounds= parseBounds(bounds_f);

    float*** tetrahedrons_list= new float**[N_ATOMS];
    int i_wat= 0;
    for(int a= 0; a < N_ATOMS; a++) {
        if(!study_molecules[a]) continue;
        Vector center= getPos(atom_list[a]);
        std::vector<Vector> sites= getPerfectTetrahedron(center,getPos(atom_list[a+1]),getPos(atom_list[a+2]),bounds);
        
        tetrahedrons_list[i_wat]= new float*[5];

        tetrahedrons_list[i_wat][0]= new float[3];
        tetrahedrons_list[i_wat][0][0]= center.x;
        tetrahedrons_list[i_wat][0][1]= center.y;
        tetrahedrons_list[i_wat][0][2]= center.z;

        for(std::size_t i_site= 0; i_site < sites.size(); i_site++) {
            tetrahedrons_list[i_wat][i_site+1]= new float[3];
            tetrahedrons_list[i_wat][i_site+1][0]= sites[i_site].x;
            tetrahedrons_list[i_wat][i_site+1][1]= sites[i_site].y;
            tetrahedrons_list[i_wat][i_site+1][2]= sites[i_site].z;
        }
        i_wat++;
    }
    return tetrahedrons_list;
}

extern "C" float* V4S(Atom* atom_list, int N_ATOMS, float* bounds_f, bool* study_molecules, const char* type_O, const int atoms_per_water_molecule, const float R_CUT_OFF) {
    Vector bounds= parseBounds(bounds_f);

    float* V4S_list= new float[N_ATOMS];
    int idx_water= 0;
    for(int a= 0; a < N_ATOMS; a++) {
        if(!study_molecules[a]) continue;
        std::vector<Vector> sites= getPerfectTetrahedron(getPos(atom_list[a]),getPos(atom_list[a+1]),getPos(atom_list[a+2]),bounds);
        std::vector<float> sum_site= getInteractionPerSite(a,atom_list,sites,N_ATOMS,atoms_per_water_molecule,type_O,bounds,R_CUT_OFF);
        
        std::sort(sum_site.data(),sum_site.data()+sites.size());
        V4S_list[idx_water++]= sum_site[3];
    }
    return V4S_list;
}

extern "C" float*** ViSPoints(Atom* atom_list, int N_ATOMS, float* bounds_f, bool* study_molecules, const char* type_O, const int atoms_per_water_molecule, const float R_CUT_OFF) {
    Vector bounds= parseBounds(bounds_f);

    float*** V4S_list= new float**[N_ATOMS];
    int idx_water= 0;
    for(int a= 0; a < N_ATOMS; a++) {
        if(!study_molecules[a]) continue;
        std::vector<Vector> sites= getPerfectTetrahedron(getPos(atom_list[a]),getPos(atom_list[a+1]),getPos(atom_list[a+2]),bounds);
        std::vector<float> sum_site= getInteractionPerSite(a,atom_list,sites,N_ATOMS,atoms_per_water_molecule,type_O,bounds,R_CUT_OFF);
        
        float xs[4]; float ys[4]; float zs[4];
        for(std::size_t i_site= 0; i_site < sites.size(); i_site++) {
            xs[i_site]= sites[i_site].x;
            ys[i_site]= sites[i_site].y;
            zs[i_site]= sites[i_site].z;
        }
        float** data_tetrahedron= cosort(sum_site,xs,ys,zs);
        V4S_list[idx_water++]= data_tetrahedron;
    }
    return V4S_list;
}

extern "C" void freeMemory(float* pointer) {
    delete[] pointer;
}

extern "C" void freeMemoryTet(float*** pointer, const int N_ATOMS, const int N_POINTS) {
    for(int i= 0; i < N_ATOMS; i++) {
        for(int j= 0; j < N_POINTS; j++)
            delete[] pointer[i][j];
        delete[] pointer[i];
    }
    delete[] pointer;
}
