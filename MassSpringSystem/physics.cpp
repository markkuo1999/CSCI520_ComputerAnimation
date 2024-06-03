/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "physics.h"
#include <iostream>
#include <cmath>
using namespace std;

struct ForceFieldPointOnAxis* firstpoint = NULL;
bool flag = true;
double forcefieldmax;

struct ForceFieldPointOnAxis
{
	double p;
	struct ForceFieldPointOnAxis* next;
};

/* insert new node at the end of the linked list */
/* represent all the points of the forcefield on one axis */
void append(struct ForceFieldPointOnAxis** head, double node_data)
{
	/* 1. create and allocate node */
	struct ForceFieldPointOnAxis* newNode = new ForceFieldPointOnAxis;

	struct ForceFieldPointOnAxis* last = *head; /* used in step 5*/

	/* 2. assign data to the node */
	newNode->p = node_data;

	/* 3. set next pointer of new node to null as its the last node*/
	newNode->next = NULL;

	/* 4. if list is empty, new node becomes first node */
	if (*head == NULL)
	{
		*head = newNode;
		return;
	}

	/* 5. Else traverse till the last node */
	while (last->next != NULL)
		last = last->next;

	/* 6. Change the next of last node */
	last->next = newNode;
	return;
}

double CountFloor(double c, struct world* jello, struct ForceFieldPointOnAxis* firstpoint)
{
	struct ForceFieldPointOnAxis* tmp = firstpoint;
	while (tmp->next != NULL)
	{
		if (c > tmp->p && c < tmp->next->p)
		{
			return tmp->p;
		}
		tmp = tmp->next;
	}
}

double CountCeil(double c, struct world* jello, struct ForceFieldPointOnAxis* firstpoint)
{
	struct ForceFieldPointOnAxis* tmp = firstpoint;
	while (tmp->next != NULL)
	{
		if (c > tmp->p && c < tmp->next->p)
		{
			return tmp->next->p;
		}
		tmp = tmp->next;
	}
}

struct point ComputeForce(int ia, int ja, int ka, int ib, int jb, int kb, struct point a, struct world * jello)
{
	double F_Hook[3], F_Damp[3], F_ForceField[3], L[3], L_length, Rest_Length;

	F_ForceField[0] = F_ForceField[1] = F_ForceField[2] = 0;
	F_Hook[0] = F_Hook[1] = F_Hook[2] = 0;
	F_Damp[0] = F_Damp[1] = F_Damp[2] = 0;
	L[0] = jello->p[ia][ja][ka].x - jello->p[ib][jb][kb].x;
	L[1] = jello->p[ia][ja][ka].y - jello->p[ib][jb][kb].y;
	L[2] = jello->p[ia][ja][ka].z - jello->p[ib][jb][kb].z;
	L_length = sqrt(L[0] * L[0] + L[1] * L[1] + L[2] * L[2]);
	Rest_Length = sqrt((ia-ib)*(ia-ib) + (ja-jb)*(ja-jb) + (ka-kb)*(ka-kb))/7;
	F_Hook[0] = -jello->kElastic * (L_length - Rest_Length) * (L[0] / L_length);
	F_Hook[1] = -jello->kElastic * (L_length - Rest_Length) * (L[1] / L_length);
	F_Hook[2] = -jello->kElastic * (L_length - Rest_Length) * (L[2] / L_length);
	F_Damp[0] = -jello->dElastic * (((jello->v[ia][ja][ka].x - jello->v[ib][jb][kb].x) * L[0]) / L_length) * (L[0] / L_length);
	F_Damp[1] = -jello->dElastic * (((jello->v[ia][ja][ka].y - jello->v[ib][jb][kb].y) * L[1]) / L_length) * (L[1] / L_length);
	F_Damp[2] = -jello->dElastic * (((jello->v[ia][ja][ka].z - jello->v[ib][jb][kb].z) * L[2]) / L_length) * (L[2] / L_length);
	a.x = 1*(F_Hook[0] + F_Damp[0] + F_ForceField[0]) / jello->mass;
	a.y = 1*(F_Hook[1] + F_Damp[1] + F_ForceField[1]) / jello->mass;
	a.z = 1*(F_Hook[2] + F_Damp[2] + F_ForceField[2]) / jello->mass;
	//cout << a.x << " " << a.y << " " << a.z << endl;
	return a;
}

struct point ComputeCollisionForce(int ia, int ja, int ka, double ib, double jb, double kb, double Distance, struct point a, struct world* jello)
{
	double F_Hook[3], F_Damp[3], F_ForceField[3], L[3], L_length, Rest_Length;
	
	F_ForceField[0] = F_ForceField[1] = F_ForceField[2] = 0;
	F_Hook[0] = F_Hook[1] = F_Hook[2] = 0;
	F_Damp[0] = F_Damp[1] = F_Damp[2] = 0;
	L[0] = jello->p[ia][ja][ka].x - ib;
	L[1] = jello->p[ia][ja][ka].y - jb;
	L[2] = jello->p[ia][ja][ka].z - kb;
	//cout << L[0] << endl;
	L_length = sqrt(L[0] * L[0] + L[1] * L[1] + L[2] * L[2]);
	Rest_Length = 0;
	F_Hook[0] = -jello->kCollision * (L_length - Rest_Length) * (L[0] / L_length);
	F_Hook[1] = -jello->kCollision * (L_length - Rest_Length) * (L[1] / L_length);
	F_Hook[2] = -jello->kCollision * (L_length - Rest_Length) * (L[2] / L_length);
	//cout << F_Hook[0] << " " << F_Hook[1] << " " << F_Hook[2] << endl;
	F_Damp[0] = -jello->dCollision * (((jello->v[ia][ja][ka].x - 0) * L[0]) / L_length) * (L[0] / L_length);
	F_Damp[1] = -jello->dCollision * (((jello->v[ia][ja][ka].y - 0) * L[1]) / L_length) * (L[1] / L_length);
	F_Damp[2] = -jello->dCollision * (((jello->v[ia][ja][ka].z - 0) * L[2]) / L_length) * (L[2] / L_length);
	//F_Damp[0] = -jello->dCollision * jello->v[ia][ja][ka].x;
	//F_Damp[1] = -jello->dCollision * jello->v[ia][ja][ka].y;
	//F_Damp[2] = -jello->dCollision * jello->v[ia][ja][ka].z;
	//cout << F_Damp[0] << " " << F_Damp[1] << " " << F_Damp[2] << endl;
	a.x = (F_Hook[0] + F_Damp[0] + F_ForceField[0]) / jello->mass;
	a.y = (F_Hook[1] + F_Damp[1] + F_ForceField[1]) / jello->mass;
	a.z = (F_Hook[2] + F_Damp[2] + F_ForceField[2]) / jello->mass;
	//cout << a.x << " " << a.y << " " << a.z << endl;
	return a;
}

struct ForcePoint
{
	double x;
	double y;
	double z;
};

/* Computes acceleration to every control point of the jello cube, 
   which is in state given by 'jello'.
   Returns result in array 'a'. */
void computeAcceleration(struct world * jello, struct point a[8][8][8])
{
  /* for you to implement ... */
  int i, j, k;
  int b, c, d;
  double Distance;
  double alpha, beta, gama;
  struct ForcePoint forcepoint[8], TotalForceField;

  for (i = 0; i < 8; i++)
  {
	  forcepoint[i].x = 0;
	  forcepoint[i].y = 0;
	  forcepoint[i].z = 0;
  }

  for (i = 0; i <= 7; i++)
	  for (j = 0; j <= 7; j++)
		  for (k = 0; k <= 7; k++)
		  {	
			  a[i][j][k].x = a[i][j][k].y = a[i][j][k].z = 0;
			  /*Structral Springs*/
			  if (i + 1 <= 7)
			  {
				  a[i][j][k].x += ComputeForce(i, j, k, i + 1, j, k, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeForce(i, j, k, i + 1, j, k, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeForce(i, j, k, i + 1, j, k, a[i][j][k], jello).z;
			  }
			  if (i - 1 >= 0)
			  {
				  a[i][j][k].x += ComputeForce(i, j, k, i - 1, j, k, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeForce(i, j, k, i - 1, j, k, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeForce(i, j, k, i - 1, j, k, a[i][j][k], jello).z;
			  } 
			  if (j + 1 <= 7)
			  {
				  a[i][j][k].x += ComputeForce(i, j, k, i, j + 1, k, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeForce(i, j, k, i, j + 1, k, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeForce(i, j, k, i, j + 1, k, a[i][j][k], jello).z;
			  }
			  if (j - 1 >= 0)
			  {
				  a[i][j][k].x += ComputeForce(i, j, k, i, j - 1, k, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeForce(i, j, k, i, j - 1, k, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeForce(i, j, k, i, j - 1, k, a[i][j][k], jello).z;
			  }
			  if (k + 1 <= 7)
			  {
				  a[i][j][k].x += ComputeForce(i, j, k, i, j, k + 1, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeForce(i, j, k, i, j, k + 1, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeForce(i, j, k, i, j, k + 1, a[i][j][k], jello).z;
			  }
			  if (k - 1 >= 0)
			  {
				  a[i][j][k].x += ComputeForce(i, j, k, i, j, k - 1, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeForce(i, j, k, i, j, k - 1, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeForce(i, j, k, i, j, k - 1, a[i][j][k], jello).z;
			  }
			  /*Shear Springs*/
			  if (i + 1 <= 7 && j + 1 <= 7)
			  {
				  a[i][j][k].x += ComputeForce(i, j, k, i + 1, j + 1, k, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeForce(i, j, k, i + 1, j + 1, k, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeForce(i, j, k, i + 1, j + 1, k, a[i][j][k], jello).z;
			  }
			  if (i - 1 >= 0 && j + 1 <= 7)
			  {
				  a[i][j][k].x += ComputeForce(i, j, k, i - 1, j + 1, k, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeForce(i, j, k, i - 1, j + 1, k, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeForce(i, j, k, i - 1, j + 1, k, a[i][j][k], jello).z;
			  }
			  if (i - 1 >= 0 && j - 1 >= 0)
			  {
				  a[i][j][k].x += ComputeForce(i, j, k, i - 1, j - 1, k, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeForce(i, j, k, i - 1, j - 1, k, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeForce(i, j, k, i - 1, j - 1, k, a[i][j][k], jello).z;
			  }
			  if (i + 1 <= 7 && j - 1 >= 0)
			  {
				  a[i][j][k].x += ComputeForce(i, j, k, i + 1, j - 1, k, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeForce(i, j, k, i + 1, j - 1, k, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeForce(i, j, k, i + 1, j - 1, k, a[i][j][k], jello).z;
			  }
			  if (j + 1 <= 7 && k + 1 <= 7)
			  {
				  a[i][j][k].x += ComputeForce(i, j, k, i, j + 1, k + 1, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeForce(i, j, k, i, j + 1, k + 1, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeForce(i, j, k, i, j + 1, k + 1, a[i][j][k], jello).z;
			  }
			  if (j - 1 >= 0 && k + 1 <= 7)
			  {
				  a[i][j][k].x += ComputeForce(i, j, k, i, j - 1, k + 1, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeForce(i, j, k, i, j - 1, k + 1, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeForce(i, j, k, i, j - 1, k + 1, a[i][j][k], jello).z;
			  }
			  if (j - 1 >= 0 && k - 1 >= 0)
			  {
				  a[i][j][k].x += ComputeForce(i, j, k, i, j - 1, k - 1, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeForce(i, j, k, i, j - 1, k - 1, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeForce(i, j, k, i, j - 1, k - 1, a[i][j][k], jello).z;
			  }
			  if (j + 1 <= 7 && k - 1 >= 0)
			  {
				  a[i][j][k].x += ComputeForce(i, j, k, i, j + 1, k - 1, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeForce(i, j, k, i, j + 1, k - 1, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeForce(i, j, k, i, j + 1, k - 1, a[i][j][k], jello).z;
			  }
			  if (i + 1 <= 7 && k + 1 <= 7)
			  {
				  a[i][j][k].x += ComputeForce(i, j, k, i + 1, j, k + 1, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeForce(i, j, k, i + 1, j, k + 1, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeForce(i, j, k, i + 1, j, k + 1, a[i][j][k], jello).z;
			  }
			  if (i - 1 >= 0 && k + 1 <= 7)
			  {
				  a[i][j][k].x += ComputeForce(i, j, k, i - 1, j, k + 1, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeForce(i, j, k, i - 1, j, k + 1, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeForce(i, j, k, i - 1, j, k + 1, a[i][j][k], jello).z;
			  }
			  if (i - 1 >= 0 && k - 1 >= 0)
			  {
				  a[i][j][k].x += ComputeForce(i, j, k, i - 1, j, k - 1, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeForce(i, j, k, i - 1, j, k - 1, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeForce(i, j, k, i - 1, j, k - 1, a[i][j][k], jello).z;
			  }
			  if (i + 1 <= 7 && k - 1 >= 0)
			  {
				  a[i][j][k].x += ComputeForce(i, j, k, i + 1, j, k - 1, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeForce(i, j, k, i + 1, j, k - 1, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeForce(i, j, k, i + 1, j, k - 1, a[i][j][k], jello).z;
			  }
			  if (i + 1 <= 7 && j + 1 <= 7 && k + 1 <= 7)
			  {
				  a[i][j][k].x += ComputeForce(i, j, k, i + 1, j + 1, k + 1, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeForce(i, j, k, i + 1, j + 1, k + 1, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeForce(i, j, k, i + 1, j + 1, k + 1, a[i][j][k], jello).z;
			  }
			  if (i - 1 >= 0 && j + 1 <= 7 && k + 1 <= 7)
			  {
				  a[i][j][k].x += ComputeForce(i, j, k, i - 1, j + 1, k + 1, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeForce(i, j, k, i - 1, j + 1, k + 1, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeForce(i, j, k, i - 1, j + 1, k + 1, a[i][j][k], jello).z;
			  }
			  if (i - 1 >= 0 && j - 1 >= 0 && k + 1 <= 7)
			  {
				  a[i][j][k].x += ComputeForce(i, j, k, i - 1, j - 1, k + 1, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeForce(i, j, k, i - 1, j - 1, k + 1, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeForce(i, j, k, i - 1, j - 1, k + 1, a[i][j][k], jello).z;
			  }
			  if (i + 1 <= 7 && j - 1 >= 0 && k + 1 <= 7)
			  {
				  a[i][j][k].x += ComputeForce(i, j, k, i + 1, j - 1, k + 1, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeForce(i, j, k, i + 1, j - 1, k + 1, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeForce(i, j, k, i + 1, j - 1, k + 1, a[i][j][k], jello).z;
			  }
			  if (i + 1 <= 7 && j + 1 <= 7 && k - 1 >= 0)
			  {
				  a[i][j][k].x += ComputeForce(i, j, k, i + 1, j + 1, k - 1, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeForce(i, j, k, i + 1, j + 1, k - 1, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeForce(i, j, k, i + 1, j + 1, k - 1, a[i][j][k], jello).z;
			  }
			  if (i - 1 >=0 && j + 1 <= 7 && k - 1 >= 0)
			  {
				  a[i][j][k].x += ComputeForce(i, j, k, i - 1, j + 1, k - 1, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeForce(i, j, k, i - 1, j + 1, k - 1, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeForce(i, j, k, i - 1, j + 1, k - 1, a[i][j][k], jello).z;
			  }
			  if (i - 1 >= 0 && j - 1 >= 0 && k - 1 >= 0)
			  {
				  a[i][j][k].x += ComputeForce(i, j, k, i - 1, j - 1, k - 1, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeForce(i, j, k, i - 1, j - 1, k - 1, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeForce(i, j, k, i - 1, j - 1, k - 1, a[i][j][k], jello).z;
			  }
			  if (i + 1 <= 7 && j - 1 >= 0 && k - 1 >= 0)
			  {
				  a[i][j][k].x += ComputeForce(i, j, k, i + 1, j - 1, k - 1, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeForce(i, j, k, i + 1, j - 1, k - 1, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeForce(i, j, k, i + 1, j - 1, k - 1, a[i][j][k], jello).z;
			  }
			  /*Bend Springs*/
			  if (i + 2 <= 7)
			  {
				  a[i][j][k].x += ComputeForce(i, j, k, i + 2, j, k, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeForce(i, j, k, i + 2, j, k, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeForce(i, j, k, i + 2, j, k, a[i][j][k], jello).z;
			  }
			  if (i - 2 >= 0)
			  {
				  a[i][j][k].x += ComputeForce(i, j, k, i - 2, j, k, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeForce(i, j, k, i - 2, j, k, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeForce(i, j, k, i - 2, j, k, a[i][j][k], jello).z;
			  }
			  if (j + 2 <= 7)
			  {
				  a[i][j][k].x += ComputeForce(i, j, k, i, j + 2, k, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeForce(i, j, k, i, j + 2, k, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeForce(i, j, k, i, j + 2, k, a[i][j][k], jello).z;
			  }
			  if (j - 2 >= 0)
			  {
				  a[i][j][k].x += ComputeForce(i, j, k, i, j - 2, k, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeForce(i, j, k, i, j - 2, k, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeForce(i, j, k, i, j - 2, k, a[i][j][k], jello).z;
			  }
			  if (k + 2 <= 7)
			  {
				  a[i][j][k].x += ComputeForce(i, j, k, i, j, k + 2, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeForce(i, j, k, i, j, k + 2, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeForce(i, j, k, i, j, k + 2, a[i][j][k], jello).z;
			  }
			  if (k - 2 >= 0)
			  {
				  a[i][j][k].x += ComputeForce(i, j, k, i, j, k - 2, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeForce(i, j, k, i, j, k - 2, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeForce(i, j, k, i, j, k - 2, a[i][j][k], jello).z;
			  }
			  //cout << 4 / double(jello->resolution) << endl;
			  double length_of_grid = 4 / (double(jello->resolution)-1);
			  //cout << length_of_grid << endl;
			  //cout << CountFloor(jello->p[i][j][k].x, jello, firstpoint) << endl;
			  //cout << CountCeil(jello->p[i][j][k].x, jello, firstpoint) << endl;
			  if(jello->p[i][j][k].x > -2 && jello->p[i][j][k].x < 2 && jello->p[i][j][k].y > -2 && jello->p[i][j][k].y < 2 && jello->p[i][j][k].z > -2 && jello->p[i][j][k].z < 2)
				  if ((CountCeil(jello->p[i][j][k].x, jello, firstpoint) - CountFloor(jello->p[i][j][k].x, jello, firstpoint)) != 0 && (CountCeil(jello->p[i][j][k].y, jello, firstpoint) - CountFloor(jello->p[i][j][k].y, jello, firstpoint)) != 0 && 
					  (CountCeil(jello->p[i][j][k].z, jello, firstpoint) - CountFloor(jello->p[i][j][k].z, jello, firstpoint)) != 0)
				  {
					 
					  //cout << "ceil " << CountCeil(jello->p[i][j][k].x, jello, firstpoint) << endl;
					  //cout << "floor " << CountFloor(jello->p[i][j][k].x, jello, firstpoint) << endl;
					  alpha = (jello->p[i][j][k].x - CountFloor(jello->p[i][j][k].x, jello, firstpoint)) / (CountCeil(jello->p[i][j][k].x, jello, firstpoint) - CountFloor(jello->p[i][j][k].x, jello, firstpoint));
					  beta = (jello->p[i][j][k].y - CountFloor(jello->p[i][j][k].y, jello, firstpoint)) / (CountCeil(jello->p[i][j][k].y, jello, firstpoint) - CountFloor(jello->p[i][j][k].y, jello, firstpoint));
					  gama = (jello->p[i][j][k].z - CountFloor(jello->p[i][j][k].z, jello, firstpoint)) / (CountCeil(jello->p[i][j][k].z, jello, firstpoint) - CountFloor(jello->p[i][j][k].z, jello, firstpoint));
					  //cout << jello->p[i][j][k].x - CountFloor(jello->p[i][j][k].x, jello, firstpoint) << endl;
					  //cout << alpha << " " << beta << " " << gama << endl;
					  
					  if (int(floor((jello->p[i][j][k].x + 2) / length_of_grid)) * jello->resolution * jello->resolution
						  + int(floor((jello->p[i][j][k].y + 2) / length_of_grid)) * jello->resolution + int(floor((jello->p[i][j][k].z + 2) / length_of_grid)) < 0)
						  cout << int(floor((jello->p[i][j][k].x + 2) / length_of_grid)) * jello->resolution * jello->resolution
						  + int(floor((jello->p[i][j][k].y + 2) / length_of_grid)) * jello->resolution + int(floor((jello->p[i][j][k].z + 2) / length_of_grid)) << endl;
					  /*A000*/
					  //cout << int(floor((jello->p[i][j][k].x + 2) / length_of_grid)) * jello->resolution * jello->resolution
						  //+ int(floor((jello->p[i][j][k].y + 2) / length_of_grid)) * jello->resolution + int(floor((jello->p[i][j][k].z + 2) / length_of_grid)) << endl;
					  forcepoint[0].x = jello->forceField[int(floor((jello->p[i][j][k].x + 2) / length_of_grid)) * jello->resolution * jello->resolution
						  + int(floor((jello->p[i][j][k].y + 2) / length_of_grid)) * jello->resolution + int(floor((jello->p[i][j][k].z + 2) / length_of_grid))].x;
					  forcepoint[0].y = jello->forceField[int(floor((jello->p[i][j][k].x + 2) / length_of_grid)) * jello->resolution * jello->resolution
						  + int(floor((jello->p[i][j][k].y + 2) / length_of_grid)) * jello->resolution + int(floor((jello->p[i][j][k].z + 2) / length_of_grid))].y;
					  forcepoint[0].z = jello->forceField[int(floor((jello->p[i][j][k].x + 2) / length_of_grid)) * jello->resolution * jello->resolution
						  + int(floor((jello->p[i][j][k].y + 2) / length_of_grid)) * jello->resolution + int(floor((jello->p[i][j][k].z + 2) / length_of_grid))].z;

					  /*A100*/
					  forcepoint[1].x = jello->forceField[int(ceil((jello->p[i][j][k].x + 2) / length_of_grid)) * jello->resolution * jello->resolution
						  + int(floor((jello->p[i][j][k].y + 2) / length_of_grid)) * jello->resolution + int(floor((jello->p[i][j][k].z + 2) / length_of_grid))].x;
					  forcepoint[1].y = jello->forceField[int(ceil((jello->p[i][j][k].x + 2) / length_of_grid)) * jello->resolution * jello->resolution
						  + int(floor((jello->p[i][j][k].y + 2) / length_of_grid)) * jello->resolution + int(floor((jello->p[i][j][k].z + 2) / length_of_grid))].y;
					  forcepoint[1].z = jello->forceField[int(ceil((jello->p[i][j][k].x + 2) / length_of_grid)) * jello->resolution * jello->resolution
						  + int(floor((jello->p[i][j][k].y + 2) / length_of_grid)) * jello->resolution + int(floor((jello->p[i][j][k].z + 2) / length_of_grid))].z;



					  /*A110*/
					  forcepoint[2].x = jello->forceField[int(ceil((jello->p[i][j][k].x + 2) / length_of_grid)) * jello->resolution * jello->resolution
						  + int(ceil((jello->p[i][j][k].y + 2) / length_of_grid)) * jello->resolution + int(floor((jello->p[i][j][k].z + 2) / length_of_grid))].x;
					  forcepoint[2].y = jello->forceField[int(ceil((jello->p[i][j][k].x + 2) / length_of_grid)) * jello->resolution * jello->resolution
						  + int(ceil((jello->p[i][j][k].y + 2) / length_of_grid)) * jello->resolution + int(floor((jello->p[i][j][k].z + 2) / length_of_grid))].y;
					  forcepoint[2].z = jello->forceField[int(ceil((jello->p[i][j][k].x + 2) / length_of_grid)) * jello->resolution * jello->resolution
						  + int(ceil((jello->p[i][j][k].y + 2) / length_of_grid)) * jello->resolution + int(floor((jello->p[i][j][k].z + 2) / length_of_grid))].z;

					  /*A010*/
					  forcepoint[3].x = jello->forceField[int(floor((jello->p[i][j][k].x + 2) / length_of_grid)) * jello->resolution * jello->resolution
						  + int(ceil((jello->p[i][j][k].y + 2) / length_of_grid)) * jello->resolution + int(floor((jello->p[i][j][k].z + 2) / length_of_grid))].x;
					  forcepoint[3].y = jello->forceField[int(floor((jello->p[i][j][k].x + 2) / length_of_grid)) * jello->resolution * jello->resolution
						  + int(ceil((jello->p[i][j][k].y + 2) / length_of_grid)) * jello->resolution + int(floor((jello->p[i][j][k].z + 2) / length_of_grid))].y;
					  forcepoint[3].z = jello->forceField[int(floor((jello->p[i][j][k].x + 2) / length_of_grid)) * jello->resolution * jello->resolution
						  + int(ceil((jello->p[i][j][k].y + 2) / length_of_grid)) * jello->resolution + int(floor((jello->p[i][j][k].z + 2) / length_of_grid))].z;

					  /*A001*/
					  forcepoint[jello->resolution/2].x = jello->forceField[int(floor((jello->p[i][j][k].x + 2) / length_of_grid)) * jello->resolution * jello->resolution
						  + int(floor((jello->p[i][j][k].y + 2) / length_of_grid)) * jello->resolution + int(ceil((jello->p[i][j][k].z + 2) / length_of_grid))].x;
					  forcepoint[jello->resolution/2].y = jello->forceField[int(floor((jello->p[i][j][k].x + 2) / length_of_grid)) * jello->resolution * jello->resolution
						  + int(floor((jello->p[i][j][k].y + 2) / length_of_grid)) * jello->resolution + int(ceil((jello->p[i][j][k].z + 2) / length_of_grid))].y;
					  forcepoint[jello->resolution/2].z = jello->forceField[int(floor((jello->p[i][j][k].x + 2) / length_of_grid)) * jello->resolution * jello->resolution
						  + int(floor((jello->p[i][j][k].y + 2) / length_of_grid)) * jello->resolution + int(ceil((jello->p[i][j][k].z + 2) / length_of_grid))].z;

					  /*A101*/
					  forcepoint[5].x = jello->forceField[int(ceil((jello->p[i][j][k].x + 2) / length_of_grid)) * jello->resolution * jello->resolution
						  + int(floor((jello->p[i][j][k].y + 2) / length_of_grid)) * jello->resolution + int(ceil((jello->p[i][j][k].z + 2) / length_of_grid))].x;
					  forcepoint[5].y = jello->forceField[int(ceil((jello->p[i][j][k].x + 2) / length_of_grid)) * jello->resolution * jello->resolution
						  + int(floor((jello->p[i][j][k].y + 2) / length_of_grid)) * jello->resolution + int(ceil((jello->p[i][j][k].z + 2) / length_of_grid))].y;
					  forcepoint[5].z = jello->forceField[int(ceil((jello->p[i][j][k].x + 2) / length_of_grid)) * jello->resolution * jello->resolution
						  + int(floor((jello->p[i][j][k].y + 2) / length_of_grid)) * jello->resolution + int(ceil((jello->p[i][j][k].z + 2) / length_of_grid))].z;

					  /*A111*/
					  forcepoint[6].x = jello->forceField[int(ceil((jello->p[i][j][k].x + 2) / length_of_grid)) * jello->resolution * jello->resolution
						  + int(ceil((jello->p[i][j][k].y + 2) / length_of_grid)) * jello->resolution + int(ceil((jello->p[i][j][k].z + 2) / length_of_grid))].x;
					  forcepoint[6].y = jello->forceField[int(ceil((jello->p[i][j][k].x + 2) / length_of_grid)) * jello->resolution * jello->resolution
						  + int(ceil((jello->p[i][j][k].y + 2) / length_of_grid)) * jello->resolution + int(ceil((jello->p[i][j][k].z + 2) / length_of_grid))].y;
					  forcepoint[6].z = jello->forceField[int(ceil((jello->p[i][j][k].x + 2) / length_of_grid)) * jello->resolution * jello->resolution
						  + int(ceil((jello->p[i][j][k].y + 2) / length_of_grid)) * jello->resolution + int(ceil((jello->p[i][j][k].z + 2) / length_of_grid))].z;

					  /*A011*/
					  forcepoint[7].x = jello->forceField[int(floor((jello->p[i][j][k].x + 2) / length_of_grid)) * jello->resolution * jello->resolution
						  + int(ceil((jello->p[i][j][k].y + 2) / length_of_grid)) * jello->resolution + int(ceil((jello->p[i][j][k].z + 2) / length_of_grid))].x;
					  forcepoint[7].y = jello->forceField[int(floor((jello->p[i][j][k].x + 2) / length_of_grid)) * jello->resolution * jello->resolution
						  + int(ceil((jello->p[i][j][k].y + 2) / length_of_grid)) * jello->resolution + int(ceil((jello->p[i][j][k].z + 2) / length_of_grid))].y;
					  forcepoint[7].z = jello->forceField[int(floor((jello->p[i][j][k].x + 2) / length_of_grid)) * jello->resolution * jello->resolution
						  + int(ceil((jello->p[i][j][k].y + 2) / length_of_grid)) * jello->resolution + int(ceil((jello->p[i][j][k].z + 2) / length_of_grid))].z;
		
					  TotalForceField.x = (1 - alpha) * (1 - beta) * (1 - gama) * forcepoint[0].x + (alpha) * (1 - beta) * (1 - gama) * forcepoint[1].x + (alpha) * (beta) * (1 - gama) * forcepoint[2].x
						  + (1 - alpha) * (beta) * (1 - gama) * forcepoint[3].x + (1 - alpha) * (1 - beta) * (gama)*forcepoint[4].x + (alpha) * (1 - beta) * (gama)*forcepoint[5].x
						  + (alpha) * (beta) * (gama)*forcepoint[6].x + (1 - alpha) * (beta) * (gama)*forcepoint[7].x;
					  TotalForceField.y = (1 - alpha) * (1 - beta) * (1 - gama) * forcepoint[0].y + (alpha) * (1 - beta) * (1 - gama) * forcepoint[1].y + (alpha) * (beta) * (1 - gama) * forcepoint[2].y
						  + (1 - alpha) * (beta) * (1 - gama) * forcepoint[3].y + (1 - alpha) * (1 - beta) * (gama)*forcepoint[4].y + (alpha) * (1 - beta) * (gama)*forcepoint[5].y
						  + (alpha) * (beta) * (gama)*forcepoint[6].y + (1 - alpha) * (beta) * (gama)*forcepoint[7].y;
					  TotalForceField.z = (1 - alpha) * (1 - beta) * (1 - gama) * forcepoint[0].z + (alpha) * (1 - beta) * (1 - gama) * forcepoint[1].z + (alpha) * (beta) * (1 - gama) * forcepoint[2].z
						  + (1 - alpha) * (beta) * (1 - gama) * forcepoint[3].z + (1 - alpha) * (1 - beta) * (gama)*forcepoint[4].z + (alpha) * (1 - beta) * (gama)*forcepoint[5].z
						  + (alpha) * (beta) * (gama)*forcepoint[6].z + (1 - alpha) * (beta) * (gama)*forcepoint[7].z;
					  //cout << TotalForceField.x << endl;

					  
					  a[i][j][k].x += 1*TotalForceField.x / jello->mass;
					  a[i][j][k].y += 1*TotalForceField.y / jello->mass;
					  a[i][j][k].z += 1*TotalForceField.z / jello->mass;
				  }
			  
				  
			  

			  /*Collision Springs*/
			  /*Left Face*/
			  Distance = 0 * jello->p[i][j][k].x + 1 * jello->p[i][j][k].y + 0 * jello->p[i][j][k].z + 2;
			  if (Distance < 0) {
				  a[i][j][k].x += ComputeCollisionForce(i, j, k, jello->p[i][j][k].x, -2, jello->p[i][j][k].z, Distance, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeCollisionForce(i, j, k, jello->p[i][j][k].x, -2, jello->p[i][j][k].z, Distance, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeCollisionForce(i, j, k, jello->p[i][j][k].x, -2, jello->p[i][j][k].z, Distance, a[i][j][k], jello).z;
			  }
			  /*Right Face*/
			  Distance = 0 * jello->p[i][j][k].x + 1 * jello->p[i][j][k].y + 0 * jello->p[i][j][k].z - 2;
			  if (Distance > 0) {
				  a[i][j][k].x += ComputeCollisionForce(i, j, k, jello->p[i][j][k].x, 2, jello->p[i][j][k].z, Distance, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeCollisionForce(i, j, k, jello->p[i][j][k].x, 2, jello->p[i][j][k].z, Distance, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeCollisionForce(i, j, k, jello->p[i][j][k].x, 2, jello->p[i][j][k].z, Distance, a[i][j][k], jello).z;
			  }
			  /*Back Face*/
			  Distance = 1 * jello->p[i][j][k].x + 0 * jello->p[i][j][k].y + 0 * jello->p[i][j][k].z + 2;
			  if (Distance < 0) {
				  a[i][j][k].x += ComputeCollisionForce(i, j, k, -2, jello->p[i][j][k].y, jello->p[i][j][k].z, Distance, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeCollisionForce(i, j, k, -2, jello->p[i][j][k].y, jello->p[i][j][k].z, Distance, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeCollisionForce(i, j, k, -2, jello->p[i][j][k].y, jello->p[i][j][k].z, Distance, a[i][j][k], jello).z;
			  }
			  /*Front Face*/
			  Distance = 1 * jello->p[i][j][k].x + 0 * jello->p[i][j][k].y + 0 * jello->p[i][j][k].z - 2;
			  if (Distance > 0) {
				  a[i][j][k].x += ComputeCollisionForce(i, j, k, 2, jello->p[i][j][k].y, jello->p[i][j][k].z, Distance, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeCollisionForce(i, j, k, 2, jello->p[i][j][k].y, jello->p[i][j][k].z, Distance, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeCollisionForce(i, j, k, 2, jello->p[i][j][k].y, jello->p[i][j][k].z, Distance, a[i][j][k], jello).z;
			  }
			  /*Top Face*/
			  Distance = 0 * jello->p[i][j][k].x + 0 * jello->p[i][j][k].y + 1 * jello->p[i][j][k].z - 2;
			  if (Distance > 0) {
				  a[i][j][k].x += ComputeCollisionForce(i, j, k, jello->p[i][j][k].x, jello->p[i][j][k].y, 2, Distance, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeCollisionForce(i, j, k, jello->p[i][j][k].x, jello->p[i][j][k].y, 2, Distance, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeCollisionForce(i, j, k, jello->p[i][j][k].x, jello->p[i][j][k].y, 2, Distance, a[i][j][k], jello).z;
			  }
			  /*Bottom Face*/
			  Distance = 0 * jello->p[i][j][k].x + 0 * jello->p[i][j][k].y + 1 * jello->p[i][j][k].z + 2;
			  if (Distance < 0) {
				  a[i][j][k].x += ComputeCollisionForce(i, j, k, jello->p[i][j][k].x, jello->p[i][j][k].y, -2, Distance, a[i][j][k], jello).x;
				  a[i][j][k].y += ComputeCollisionForce(i, j, k, jello->p[i][j][k].x, jello->p[i][j][k].y, -2, Distance, a[i][j][k], jello).y;
				  a[i][j][k].z += ComputeCollisionForce(i, j, k, jello->p[i][j][k].x, jello->p[i][j][k].y, -2, Distance, a[i][j][k], jello).z;
			  }
			  //cout << a[i][j][k].x << " " << a[i][j][k].y << " " << a[i][j][k].z << endl;
		  }
}

/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world * jello)
{
  int i,j,k;
  point a[8][8][8];

 
  double position = -2;
  if (flag == true)
	  for (int j = 0; j < jello->resolution; j++)
	  {
		  append(&firstpoint, position);
		  /* coordinate of the force field points on one axis*/
		  position += 4 / (double(jello->resolution) - 1);

	  }
  flag = false;
  computeAcceleration(jello, a);
  
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
        jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
        jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
        jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;
        jello->v[i][j][k].x += jello->dt * a[i][j][k].x;
        jello->v[i][j][k].y += jello->dt * a[i][j][k].y;
        jello->v[i][j][k].z += jello->dt * a[i][j][k].z;

      }
}

/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */
void RK4(struct world * jello)
{
  point F1p[8][8][8], F1v[8][8][8], 
        F2p[8][8][8], F2v[8][8][8],
        F3p[8][8][8], F3v[8][8][8],
        F4p[8][8][8], F4v[8][8][8];

  point a[8][8][8];


  struct world buffer;

  int i,j,k;

  buffer = *jello; // make a copy of jello

  double position = -2;
  if (flag == true)
	  for (int j = 0; j < jello->resolution; j++)
	  {
		  //cout << position << endl;
		  append(&firstpoint, position);
		  position += 4 / (double(jello->resolution) - 1);

	  }
  flag = false;
  computeAcceleration(jello, a);
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         pMULTIPLY(jello->v[i][j][k],jello->dt,F1p[i][j][k]);
         pMULTIPLY(a[i][j][k],jello->dt,F1v[i][j][k]);
         pMULTIPLY(F1p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F1v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F2p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F2p[i][j][k]);
         // F2v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F2v[i][j][k]);
         pMULTIPLY(F2p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F2v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F3p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F3v[i][j][k]);
         pMULTIPLY(F3p[i][j][k],1.0,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],1.0,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }
         
  computeAcceleration(&buffer, a);


  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F4p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F4v[i][j][k]);

         pMULTIPLY(F2p[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3p[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1p[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4p[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->p[i][j][k],jello->p[i][j][k]);

         pMULTIPLY(F2v[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4v[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->v[i][j][k],jello->v[i][j][k]);
      }

  return;  
}
