/*

In-Place Loop Subdivision
Submission Author: John Varghese George
Description: 
This program subdivides a given input mesh using Loop Subdivision scheme
and writes the output to the given output file. The program subdivides
dynamically by using edge split and edge flip operations on the orginal
mesh object.

*/
#include <fstream>
#include <vector>

#include "OBJFileReader.h"
#include "SolidDelegate.h"
#include "Solid.h"
#include "iterators.h"
#include <list>

using namespace MeshLib;

int main(int argc, char *argv[])
{
	// Read in the obj file
	Solid mesh;
	OBJFileReader of;
	std::ifstream in(argv[1]);
	of.readToSolid(&mesh, in);

	SolidDelegate delegate;

	SolidVertexIterator viter(&mesh);
	std::map<Vertex*,Point*> evenMap;
	int max_vid=-1;
	Point *pt;

	// Calculate Even vertex points
	std::cout<<"Calculating Even Vertex Positions..."<<std::endl;
	for(;!viter.end();++viter){
		Vertex *v = *viter;
		pt = new Point();
		int n = 0;
		max_vid = ( max_vid > v->id() )? max_vid: v->id(); 	//Update max vertex id as we go
		if(!(v->boundary())){							   	//Interior Vertex Case
			HalfEdge *he = v->halfedge();
			do{
				Vertex *adj = he->source();					//Get neighbour
				*pt += adj->point();					  	//Add neighbour point
				++n;									  	//Keep track of number of neighbours
				he = he->clw_rotate_about_target();		  	//Go to next neighbour
			}while(he != v->halfedge());                  	//Loop until back to start
			double beta = 0.375 + 0.25*cos(2*M_PI/n);     	//Compute beta coeff for loop subdiv
			beta = (1.0/n)*(0.625-(beta*beta));
			*pt *= beta;
			*pt += (v->point())*(1.0-n*beta);             	//Add contribution of original vertex
		}else {																			//Boundary Vertex Case
			Vertex *boundary_neighbour1 = v->most_ccw_in_halfedge()->source();			//Get first connected boundary vertex
			Vertex *boundary_neighbour2 = v->most_clw_out_halfedge()->target();			//Get second connected boundary vertex
			*pt += (boundary_neighbour1->point() + boundary_neighbour2->point())*0.125;	//Add neighbour contributions
			*pt += (v->point())*0.75;													//Add original vertex contribution
		}
		evenMap[v] = pt;								  	//Store new point in hashmap for future update
	}

	SolidEdgeIterator eiter(&mesh);
	std::map<Edge*,Point*> oddMap;

	// Calculate Odd vertex points
	std::cout<<"Calculating Odd Vertex Positions..."<<std::endl;
	for(;!eiter.end();++eiter){
		Edge *e = *eiter;
		Vertex *a,*b,*c,*d;
		pt = new Point();
		e->get_vertices(a,b);                     	//Get edge points
		if(!(e->boundary())){															//Interior Vertex Case
			c = e->halfedge(0)->he_next()->target();  									//Get first corner point
			d = e->halfedge(1)->he_next()->target();									//Get second corner point
			*pt += (a->point() + b->point())*0.375 + (c->point() + d->point())*0.125;	//Compute new point using loop subdiv
		}else {
			*pt += (a->point() + b->point())*0.5;  //Boundary Vertex case
		}
		oddMap[e] = pt;			//Store new point in hashmap for future update
	}

	//Update even vertices
	std::cout<<"Updating Even Vertices..."<<std::endl;
	viter.reset();
	for(;!viter.end();++viter){
		Vertex *v = *viter;
		v->point() = *(evenMap[v]);
		delete (evenMap[v]);			//Clear memory from map since data already updated
	}
	evenMap.clear();

	std::list<Edge*> edgesToSplit;		//Create a separate edge iteratable for
	eiter.reset();						//current edges since new edges will be added when split
	for(;!eiter.end();++eiter){
		Edge *e = *eiter;
		edgesToSplit.push_back(e);
	}

	//Spilt edges to create odd vertices
	std::list<Edge*> edgesToFlip;		//Keep track of new edges that need to be corrected by flipping

	std::cout<<"Splitting edges and creating odd vertices..."<<std::endl;
	while(!edgesToSplit.empty()){
		Edge *e = edgesToSplit.front();							//Get the next edge from list
		Vertex *a,*b,*c,*d,*new_v;
		e->get_vertices(a,b);                     				//Get edge points
		if(e->halfedge(0))										//Check if halfedge exists (for boundary case)
			c = e->halfedge(0)->he_next()->target();  			//Get opposite corner
		else
			c = NULL;											//else set to NULL explicitly
		if(e->halfedge(1))										//Same check but for opposite halfedge
			d = e->halfedge(1)->he_next()->target();			//Get other corner
		else
			d = NULL;											//else set to NULL
		new_v = delegate.edgeSplit(&mesh,e);					//Perform edge split
		new_v->point() = *(oddMap[e]);  						//Update odd vertex position
		delete (oddMap[e]);										//Clear memory from map since data already updated
		if((c)&&(c->id() <= max_vid))							//Check if created edges to corners need correction
			edgesToFlip.push_back(mesh.vertexEdge(new_v,c));	//after splitting is complete.
		if((d)&&(d->id() <= max_vid))							//Add those edges to flip list
			edgesToFlip.push_back(mesh.vertexEdge(new_v,d));	//Null checks here identify new boundary edges
		edgesToSplit.pop_front();
	}
	oddMap.clear();

	//Edge flip code
	//(Could not find function to do this in library so implemented here)
	std::cout<<"Correcting select edges by flipping..."<<std::endl;
	while(!edgesToFlip.empty()){
		Edge *e = edgesToFlip.front();				//Get next edge in list
		Vertex *a,*b,*c,*d;
		e->get_vertices(a,b);						//Get edge vertices
		c = e->halfedge(0)->he_next()->target();	//Get first corner
		d = e->halfedge(1)->he_next()->target();	//Get second corner
		Face *f1,*f2;
		int f1_id,f2_id,v[3];
		f1 = e->halfedge(0)->face();				//Get first face with current edge
		f1_id = f1->id();
		f2 = e->halfedge(1)->face();				//Get second face with current edge
		f2_id = f2->id();							//Both faces exist since we already check for boundary before
		delegate.removeFace(&mesh,f1);				//Remove both faces
		delegate.removeFace(&mesh,f2);
		v[0] = a->id();								//Create new faces with flipped edge
		v[1] = d->id();								//Keep original face ids for new faces
		v[2] = c->id();
		delegate.createFace(&mesh,v,f1_id);
		v[0] = b->id();
		v[1] = c->id();
		v[2] = d->id();
		delegate.createFace(&mesh,v,f2_id);
		edgesToFlip.pop_front();					//Remove edge from processing list
	}

	// Write out the resultant obj file
	int vObjID = 1;
	std::map<int, int> vidToObjID;

	std::ofstream os(argv[2]);

	SolidVertexIterator iter(&mesh);

	for(; !iter.end(); ++iter)
	{
		Vertex *v = *iter;
		Point p = v->point();
		os << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
		vidToObjID[v->id()] = vObjID++;
	}
	os << "# " << (unsigned int)mesh.numVertices() << " vertices" << std::endl;

	float u = 0.0, v = 0.0;
	for(iter.reset(); !iter.end(); ++iter)
	{
		Vertex *vv = *iter;
		std::string key( "uv" );
		std::string s = Trait::getTraitValue (vv->string(), key );
		if( s.length() > 0 )
		{
			sscanf( s.c_str (), "%f %f", &u, &v );
		}
		os << "vt " << u << " " << v << std::endl;
	}
	os << "# " << (unsigned int)mesh.numVertices() << " texture coordinates" << std::endl;

	SolidFaceIterator fiter(&mesh);
	for(; !fiter.end(); ++fiter)
	{
		Face *f = *fiter;
		FaceVertexIterator viter(f);
		os << "f " ;
		for(; !viter.end(); ++viter)
		{
			Vertex *v = *viter;
			os << vidToObjID[v->id()] << "/" << vidToObjID[v->id()] << " ";
		}
		os << std::endl;
	}
	os.close();
	return 0;
}