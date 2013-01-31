// Tin.Cc Main Program And Misc Small Routines. 

//  Wm. Randolph Franklin,  wrf@ecse.rpi.edu, (518) 276-6077;  Fax: -6261
//  ECSE Dept., 6026 JEC, Rensselaer Polytechnic Inst, Troy NY, 12180 USA

// Time-stamp: </wrf/c/tin/Dist/tin.cc, Tue, 30 Oct 2001, 15:46:42 EST, wrf@benvolio.ecse.rpi.edu>


// INCLUDES 


#include <fcntl.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/file.h>
#include <sys/times.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>

using namespace std;


/* SUBS FUNCTION DECLARATIONS */

short int **alloc_array(int n, char *msg);
float **alloc_float_array (int n, char *msg);
void    close2(int fout, char *msg);
void    fclose2(FILE * stream, char *msg);
FILE *fdopen2(int fd, char *type, char *msg);
int flock2(int fd, int operation, char *msg);
FILE   *fopen2(char *filename, char *code, char *msg);
char   *malloc2(unsigned, char *);
int     open2(char *filename, int code, char *msg);
int     open3(char *filename, int code, int mode, char *msg);
void    perror2(char *);
void    print_array(short int **a, int n, char *msg);
void    read2(int fin, char *buf, int nchar, char *msg);
void swap(short int& a, short int& b);
void swap(int& a, int& b);
void swap(float& a, float& b);
void    write2(int fout, char *buf, int nchar, char *msg);


/* BUILTIN FUNCTION DECLARATIONS */

// extern int abs(int);
extern int     close(int);
// extern void exit(int);
extern int     fclose(FILE *);
/* extern int flock(int,int); */
extern int     fprintf();
extern off_t lseek(int fd, off_t offset, int whence);
extern int     printf();
/* extern int     read(int, char *, int); */
/* extern void    perror(char *); */
/* extern int read(int, char *, int); */
extern int     sscanf();
/* extern int     write(int, char *, int); */


// INLINE FUNCTIONS

#ifndef INLINES

/* ABS */

inline float abs(const float a)
{ return (a>=0) ? a :  -a; }

// conflicts with C 
// inline int abs(const int a)
// { return (a>=0) ? a :  -a; }


/* MAX */

inline float max(const float x, const float y)
  { return (x > y) ? x : y; }

inline int max(const int x, const int y) 
  { return (x > y) ? x : y; }

/* MIN */

inline float min(const float x, const float y)
  { return (x < y) ? x : y; }

inline int min(const int x, const int y)
  { return (x < y) ? x : y; }


/* SWAP Swap 2 variables */

inline void 
swap(int &a, int &b)
{
   int     c;
   c = a;
   a = b;
   b = c;
}

inline void 
swap(short int &a, short int &b)
{
   short int c;
   c = a;
   a = b;
   b = c;
}

inline void 
swap(float &a, float &b)
{
   float     c;
   c = a;
   a = b;
   b = c;
}

#define INLINES
#endif  




// STRUCTS AND CLASSES NEEDED FOR CLASS DEFINITIONS

struct Triangle;


// VARIABLES NEEDED FOR CLASS DEFINITIONS

extern int npoints;
extern int ntriangles;		// Current number of triangles
extern int max_allowable_num_triangles;		// Max allowable number of triangles
extern Triangle *triangles;
extern float max_err_over_whole_TIN;		// Maximum error for whole TIN. 
extern float stop_err;		// Stop splitting when max_err_over_whole_TIN is this small.
extern int worst_triangle;	// ID of worst triangle, the one with the maximum error. 
extern int tprint;              // Amount Of Printing To Do. 0 - No Routine Printing. 
extern int nsplits;		// # splits so far. 

// Statistics

extern int n_elts_pushed_onto_heap;        // Number of times elements pushed onto errs.
extern int n_elts_popped_from_heap;		// Number of times elements popped from errs.
extern int n_elts_deleted_from_heap;	// Number of times elements deleted from errs.
extern int n_swaps_while_pushing_heap;	// Number of times elements swapped pushing.
extern int n_swaps_while_popping_heap;	// Number of times elements swapped popping.
extern int max_ever_heap_size;		// Hi-water mark for heap size.
extern int n_split_tri_list_traverses;
extern int n_new_traverses;

// ROUTINE DECLARATIONS NEEDED FOR CLASS DECLARATIONS

void Test_Tri(const int t, const char *msg);
void Test_Pt(const int t, const char *msg);


// LINKED LIST CLASSES.  
// All linked lists are allocced from one global array, pairs.

class Dotpair;
class Linked_list;
extern Dotpair *dotpairs;	// Array of dotted pairs for linked lists of points
				// in triangles.			
extern Linked_list freelist;		// Free list of unallocated dotted pairs in dotpairs.
typedef int Pairindex;
const Pairindex eolist(-1);	// End of dotted pair list. 

class Dotpair			// One dotted pair, containing a point id and an index to the
{				// next dotted pair.
public:
  int     ptid;			// ID of a point, i.e., index in points array. 
  Pairindex next;		// The remaining points.
};


class Linked_list
{
public: 
  Pairindex head;		// Index in dotpairs of the 1st element.
                                // ==eolist if linked list is empty.
  Linked_list() {head = eolist; }
  Linked_list(int a) {head = a;}       // ???

  void push(int point);		// Prefix a new point to this list

  int pop();				// Pop the first point off the list.

  void check()			// Check that there is not a loop at the start.
  {
    if (head != eolist && head == dotpairs[head].next)
      perror2("Checking list found that the first pointer loops.");
  }

  int length()           // Return number of elements in list.
  {			 // Cannot easily be recursive while catching infinite loops.
    int i=0;
    Pairindex l=head;
    while (l!=eolist)
      {
	l=dotpairs[l].next;
	i++;
	if (i > npoints) perror2("Apparent infinite loop in Linked_list.length.");
      }
    return i;    
  }
};


struct Point			// One point.  A struct is 4 bytes shorter than an array.
{
  short int x,y,z;
};


class Triangle
{
public:
  int     vert[3];		// ID Numbers of the 3 verts.  vert[i] is opposite
   			        // adj[i].  The vertices around a triangle are ordered in
				// the positive direction.  Vertices are not also stored
				// in the list of points in this triangle.
  int     adj[3];		// ID Numbers of the 3 adjacent triangles, or -1 if
				// that side is exterior.
  int     heap_loc;		// Location in heap of this triangle.  -1 if not in heap.
  int     worst_pt;		// ID of the worst point.  This is its index in the
				// global points array.   
  Linked_list   l;		// Linked list of the ID numbers of the points inside this
				// triangle.
  float   maxe;		// Its error from the plane of the triangle,
				// measured vertically.  Type is not short int ob possible
				// overflow.

  Triangle& operator=(const Triangle& t) // Copy triangle but don't update heap.
    {
	if (this == &t) return *this;
	this->vert[0] = t.vert[0];
	this->vert[1] = t.vert[1];
	this->vert[2] = t.vert[2];
	this->adj[0] = t.adj[0];
	this->adj[1] = t.adj[1];
	this->adj[2] = t.adj[2];
	this->heap_loc = t.heap_loc;
//	if (a.heap_loc >= 0) 
//	    errs.elts[a.heap_loc].tri = (this - &triangles[0])/sizeof(Triangle);

	this->worst_pt = t.worst_pt;
	this->l = t.l;
	this->maxe = t.maxe;
    }
} ;



// HEAP to store pairs of (float,int) representing a triangle error and number,
// and to return and delete the triangle with the largest error, i.e., the
// triangle listed in the pair with the largest first element.  The 1st
// element of the heap is the largest.  Each element is at least as large
// as its sons.  The sons of element #i are in locs #2i and #2i+1.  To make
// the indexing simpler, the 0th element is not used.  The interface is via
// the triangle number.

// When a triangle is split, deleted, or swapped, delete its entry here.

// Each triangle has a pointer to its entry in the heap (or -1 if the
// triangle is not in the heap).  This pointer must be updated when
// the heap entry is moved.

struct Heap_elt
{
    int tri;
    float err()
    {
	return triangles[this->tri].maxe;
    }
};

ostream& operator<<(ostream &s, Heap_elt e);


int Heap_Elt_Compare(const void *a, const void *b);


void swap(Heap_elt &a, Heap_elt &b)
{
    Test_Tri(a.tri, "1st arg to heap_elt swap");
    Test_Tri(b.tri, "2nd arg to heap_elt swap");
    swap(a.tri, b.tri);
    swap(triangles[a.tri].heap_loc, triangles[b.tri].heap_loc);
}


class Heap 
{
public:
  int alloc_size;		// Allocated size, excluding 0th element.
  int heap_size;			// Current size
  Heap_elt *elts;
  
  Heap(const int size)
  {
    elts = new Heap_elt[size+1];
    if (elts==0) perror2("Cannot alloc elts in heap.");
    heap_size = 0;
    alloc_size = size;
  }
  
  void push(const int t);
  
// Pop triangle with largest err, setting global vars worst_triangle and
// max_err_over_whole_TIN.  If no triangle, set -1 and 0.0 .
  
  void pop();

// Delete entry #i of heap.

    void del(const int i);
};


// TYPEDEFS

typedef float COORDT;


// CONSTANTS

const int DEBUG0 = 1;
const float TOLERANCE = 1.0e-6;	// Numbers Smaller Than This Are Worrisome. 
const int MSWAPQ = 100;		/* Size of queue of edges to consider for swapping.
				   This should be a little bigger than you'd think
				   since edges are appended to the end while they
				   are being processed from the front. */
const int MDELETABLE = 3;	/* Size of array of deletable, null, triangles. Need
				   never be over 3. */
const int plus1[3] = {1, 2, 0};	// Plus 1, Modulo 3 
const int minus1[3] = {2, 0, 1}; // Minus 1, Modulo 3 



// VARIABLE DECLARATIONS 

extern Point *points;		/* Array of  points.  Heights  are already scaled.
				X and y are explicit  to allow for future
				generalization. */ 
extern int tinn;		// Number Of Rows, Columns Of Data 
extern int swapq[MSWAPQ][2];	/* Queue of edges to consider to swapping. Each
				   entry is given as the IDs of the 2 adjacent
				   triangles, since the edges themselves are not
				   stored explicitly.  */
extern int nswapq;		// # edges in q now. 
extern int hiswapq;		// Hi-water size for edge swap q.
extern int nswaps;		// # swaps made so far. 
extern int deletable_triangle[MDELETABLE];	/* Flagged for being null and
						   exterior, but not yet  deleted. */
extern int ndeletable;		// Size Of Contents Of Deletable_Triangle. 
extern int tsplits;		// Total Number Of Splits To Do. 
extern int tswap;		// When to swap quadrilateral diagonals.
extern int tincrfiles;		// Whether to write incremental edge and triangle files.
extern Heap *errs;   // Heap of errors of each triangle.


// ENV OPTIONS STUFF

typedef struct                 // Type to define numerical envar options 
  {				
     char   *name;		// envar 
     int     defaul;		// default.  Note 'default' is a reserved word. 
     int     minval, maxval;	// allowable range 
     int    &var;		// variable to store value into. 
  } numerical_opt_t;

const numerical_opt_t numerical_opts[] = {{"TINN", 1, 10, 1201, tinn},
					  {"TPRINT", 0, 0, 1, tprint},
					  {"TINCRFILES", 0, 0, 1, tincrfiles},
					  {"TSWAP", 1, 0, 2, tswap},
					  {"TSPLITS", 10000, 0, 2000000, tsplits}};

const int n_numerical_opts = 5;




// ROUTINE DECLARATIONS 

float   Area(const Triangle t);
double  Read_Delta_Time(),  Read_Delta_Time2();
int     Split_Triangle();
void    Add_to_Swapq(int t, int a);
void    Alloc_Dotpairs(const int n);
void    Calc_One_Eqn_and_Err(int);
void    Calc_Plane_Eqn(const int t, float *const eqn);
void    Copy_Triangle(int to, int from);
void    Cross_Product(float *const a, const float *const b, const float *const c);
void    Cross_Product(float *const a, const int *const b, const int *const c);
void    Cross_Product(int *const a, const int *const b, const int *const c);
void    Delete_Flagged_Triangles();
void    Flag_Null_Triangle_For_Deletion(int i);
void    Get_Options();
void    Init_Vars();
void    Make_First_2_Triangles();
void    Print_Help();
void    Print_Points();
void    Print_Stats();
void    Print_Time(char *msg), Print_Time2(char *msg);
void    Print_Triangle(int);
void    Proc_Swapq();
void    Read_Data();
void    Update_Adjacency(int, int);
void    Write_Edges(const char * const);
void    Write_Pgm();
void    Write_Triangles(const char *const);
void Write_Triangle_Centers(const char *const filename);


// INLINE ROUTINES

// DOT2D  - Dot Product Of 2 2-D Vectors 

inline float Dot2D(const float *const a, const float *const b)
{
   return (a)[0] * (b)[0] + (a)[1] * (b)[1];
}


// DOT3D - DOT_PRODUCT In 3-D 

inline float Dot3D(const float *const a, const float *const b)
{
   return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline float Dot3D(const short int *const a, const float *const b)
{
   return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline float Dot3D(const Point *const a, const float *const b)
{
   return a->x * b[0] + a->y * b[1] + a->z * b[2];
}


// DET2D  - Determinant Of 2 2-D Vectors 

inline float Det2D(const float *const a, const float *const b)
{
   return (a)[0] * (b)[1] - (a)[1] * (b)[0];
}


// MINEQ  - Set the first arg to the min of the 2 args, and return it.

inline float Mineq(float &a, const float b)
{
    if (b < a) a = b;
    return a;
}

// SQUARE - Convert its COORDT  arg to float and square it.

inline float Square(const COORDT k)
{
   return (float)k * (float)k;
}




// VARIABLE DEFINITIONS 

Point   *points;		// (x,y,z) of each point.
int     npoints;
Triangle *triangles;
int     ntriangles;		// Current number of triangles
int     max_allowable_num_triangles;		// Max allowable number of triangles
float   max_err_over_whole_TIN;		// Maximum error for whole TIN.
int     worst_triangle;		// ID of worst triangle, with the maximum error. 
int     tinn;			// Number of rows, columns of data 
int     nsplits;		// # splits so far. 
int     swapq[MSWAPQ][2];	/* Queue of edges to consider to swapping.  Each
				   entry is given as the IDs of the 2 adjacent
				   triangles, since the edges themselves are not
				   stored explicitly.  */
int     nswapq;			// # edges in q now. 
int     hiswapq;		// Hi-water size for edge swap q.
int     nswaps;			// # swaps made so far. 
int     deletable_triangle[MDELETABLE];		/* Flagged for being null and
						   exterior, but not yet deleted. */
int     ndeletable;		// Size Of Contents Of Deletable_Triangle. 
int     tprint;			// Amount Of Printing To Do. 0 - No Routine Printing. 
int     tswap;			// When to swap quadrilateral diagonals.
int     tincrfiles;		// Whether to write incremental edge and triangle files.
int     tsplits;		// Total Number Of Splits To Do. 

Dotpair *dotpairs;		// Array of dotted pairs for linked lists of points
				// in triangles.			
Linked_list freelist;			// Free list of unallocated dotted pairs in dotpairs.
ofstream statfile("stats"); // Each line: # splits, # triangles, max_err_over_whole_TIN, # swaps.
Heap *errs;			// Heap of errors of each triangle.
int n_elts_pushed_onto_errs_heap;		// Number of times elements pushed onto errs.
int n_elts_popped_from_heap;			// Number of times elements popped from errs.
int n_elts_deleted_from_heap;		// Number of times elements deleted from errs.
int n_elts_pushed_onto_errs_heap_swap;		// Number of times elements swapped pushing.
int n_swaps_while_popping_heap;		// Number of times elements swapped popping.
int max_ever_heap_size;			// Hi-water mark for heap size.
int n_split_tri_list_traverses;
int n_new_traverses;

float stop_err=0.0;		// Stop splitting when max_err_over_whole_TIN is this small.


// ================================================================


/**; READ_DELTA_TIME
 * Returns time in seconds since last Read_Delta_Time.  Automatically initializes
 * itself on 1st call and returns 0.
 */

double
Read_Delta_Time()
{
   static struct tms *time_buffer;
   static time_t old_time, new_time;	// in 60ths of a second 

   double  delta;

   if (time_buffer == NULL)
     {
	time_buffer = (struct tms *) malloc((unsigned) sizeof(struct tms));
	return old_time = 0;
     }
   (void) times(time_buffer);
   new_time = time_buffer->tms_utime + time_buffer->tms_stime;
   delta = (new_time - old_time) / 60.;
   old_time = new_time;
   return delta;
}

/**; PRINT_TIME Print time since last call, with a message. */

void
Print_Time(char *msg)
{
   cout << "Incr CPU Time for " << msg << " =" << Read_Delta_Time() << endl;
}

/**; READ_DELTA_TIME2
 * Like Read_Delta_Time but keeps a second incremental count.
 */

double
Read_Delta_Time2()
{
   static struct tms *time_buffer;
   static time_t old_time, new_time;	// in 60ths of a second 

   double  delta;

   if (time_buffer == NULL)
     {
	time_buffer = (struct tms *) malloc((unsigned) sizeof(struct tms));
	return old_time = 0;
     }
   (void) times(time_buffer);
   new_time = time_buffer->tms_utime + time_buffer->tms_stime;
   delta = (new_time - old_time) / 60.;
   old_time = new_time;
   return delta;
}

/**; PRINT_TIME2 Print time since last read_delta_time2 call, with a message. */

void
Print_Time2(char *msg)
{
   cout << "CPU Time for " << msg << " =" << Read_Delta_Time2() << endl;
}



// ================================================================
/**

  This version specialized for tin thus:
  - perror2 calls Print_Stats(), which is also declared here.
  - swap heap_elt decl.
  
 ALLOC_ARRAY:  Alloc an NxN short int array as N^2 contiguous entries
 ALLOC_FLOAT_ARRAY:  Alloc an NxN float array as N^2 contiguous entries
 CLOSE2:  Like close, but checks return code
 FCLOSE2
 FDOPEN2
 FLOCK2
 FOPEN2
 MALLOC2
 MAX  - for ints and floats
 MIN  - for ints and floats
 OPEN2
 OPEN3: like open2, but with mode arg also.
 PERROR2: like perror, but exits.  Also calls Print_Stats() first.
 PRINT_ARRAY
 READ2
 SWAP - for ints, shorts, and floats
 WRITE2


*/

// exit returns up thru: 15



extern void Print_Stats();

/* 
 * ALLOC_ARRAY  Alloc an NxN short int array as N^2 contiguous entries, and return
 * its address.
 */

short int **
alloc_array(int n, char *msg)
{
   short int **p;
   short int *prow;
   int     i, n2;

   n2 = (int) n *(int) n;
   p = new short int *[n];
   if (p == 0)
     {
	perror(msg);
	exit(7);
     }
   prow = new short int[n2];
   if (prow == 0)
     {
	perror(msg);
	exit(8);
     }
   for (i = 0; i < n; i++)
      p[i] = prow + i * n;
   return p;
}


/* 
 * ALLOC_FLOAT_ARRAY  Alloc an NxN float array as N^2 contiguous entries, and return
 * its address.
 */

float **
alloc_float_array(int n, char *msg)
{
   float **p;
   float  *prow;
   int     i, n2;

   n2 = (int) n *(int) n;
   p = new float *[n];
   if (p == 0)
     {
	perror(msg);
	exit(13);
     }
   prow = new float[n2];
   if (prow == 0)
     {
	perror(msg);
	exit(14);
     }
   for (i = 0; i < n; i++)
      p[i] = prow + i * n;
   return p;
}


/* CLOSE2:  Like close, but checks return code */

void
close2(int fout, char *msg)
{
   if (close(fout) < 0)
      perror(msg);
}



/* FCLOSE2  */

void
fclose2(FILE *stream, char *msg)
{
   if (fclose(stream) != 0)
      perror(msg);
}


/* FDOPEN2 */
FILE   *
fdopen2(int fd, char *type, char *msg)
{
   FILE   *stream;
   if (NULL == (stream = fdopen(fd, type)))
     {
	perror(msg);
	exit(12);
     }
   return stream;
}


/* FLOCK2 */
/** When I try to compile this with g++, i get this error:

subs.c: In function `int flock2(int, int, char *)':
subs.c:83: type `flock' does not have a constructor
*** Error code 1
make: Fatal error: Command failed for target `subs.o'

*/

int
flock2(int fd, int operation, char *msg)
{
/* extern int flock(int,int); */
   perror("Error: flock2 under construction.");

/**
    if (0 != flock (fd, operation))
      {
	  perror (msg);
	  exit (10);
      }
    return 0;
*/
}



/* FOPEN2 */

FILE   *
fopen2(char *filename, char *code, char *msg)
{
   FILE   *stream;
   if (NULL == (stream = fopen(filename, code)))
     {
	perror(msg);
	exit(4);
     }
   return stream;
}


/* MALLOC2 */
char   *
malloc2(unsigned n, char *msg)
{
   char   *p;
   if (NULL == (p = new char[n]))
     {
	perror(msg);
	exit(9);
     }
   return p;
}


/* OPEN2 */

int
open2(char *filename, int code, char *msg)
{
   int     file;
   file = open(filename, code);
   if (file == -1)
     {
	perror(msg);
	exit(2);
     }
   return file;
}

/* OPEN3 - like open2 , but with mode arg also. */

int
open3(char *filename, int code, int mode, char *msg)
{
   int     file;
   file = open(filename, code, mode);
   if (file == -1)
     {
	perror(msg);
	exit(3);
     }
   return file;
}


/* PERROR2 */

void 
perror2(char *msg)
{
  static int loop=0;

  perror(msg);
  if (++loop > 1)
    {
      printf("Infinite loop in perror2!\n");
      exit(15);
    }
  Print_Stats();
  exit(11);
}


/* PRINT_ARRAY */

void
print_array(short int **a, int n, char *msg)
{
   int     i, j;
   printf("%s, n=%d:\n", msg, n);
   for (i = 0; i < n; i++)
     {
	printf("%d: ", i);
	for (j = 0; j < n; j++)
	   printf(" %d", a[i][j]);
	printf("\n");
     }
}


/* READ2 */
void
read2(int fin, char *buf, int nchar, char *msg)
{
   if (read(fin, buf, nchar) != nchar)
     {
	perror(msg);
	exit(5);
     }
}



/* WRITE2 */

void
write2(int fout, char *buf, int nchar, char *msg)
{
   if (write(fout, buf, nchar) != nchar)
     {
	perror(msg);
	exit(6);
     }
}




// MAIN  

main()
{
  char split_label[100];
  char filename[100];
  
   cout << endl 
        << "****************************************************************" 
        << endl 
        << "[TIN - Triangular Irregular Network program, Wm. Randolph Franklin 30 Oct 01]" 
        << endl << endl;
   system("/bin/date");
   if (!statfile) perror2("Can't open stats for writing.");
   
   Get_Options();
   Init_Vars();
   Read_Data();
   Read_Delta_Time2();
   Make_First_2_Triangles();
   errs->pop();			// Get worst triangle.


// Loop, splitting the worst triangle.  Stop after either a predefined
// number of splits, or because no more triangles are splittable.  At the
// start of each iteration, the worst triangle has already been popped off
// the heap, and is stored in worst_triangle and max_err_over_whole_TIN.

   statfile << "# Statistics file for tin:" <<endl
	    << "# One line written after each triangle split." << endl
	    << "# Each line: Number of triangles split so far, current number of triangles," << endl
	    << "#  max error over whole TIN, number of edges swapped so far." << endl;

   for (nsplits = 1; nsplits <= tsplits; nsplits++)
     {
       if (Split_Triangle() != 0) break; // No triangle splittable?
       Proc_Swapq();
       errs->pop();			// Get the new worst triangle.
       if (tprint!=0)
	   cout << "After " << nsplits << " splits, ntriangles= " << ntriangles
		<< ", max_err_over_whole_TIN=" << max_err_over_whole_TIN 
		<< ", nswaps=" << nswaps << endl;
       statfile << setw(8) << nsplits << setw(8) << ntriangles << setw(9) 
		<< max_err_over_whole_TIN << setw(8) << nswaps << endl;

       if ((nsplits & 01777) == 0)
	   {
	     n_new_traverses = n_split_tri_list_traverses - n_new_traverses;
	     //	     cout << n_new_traverses << ' ';
	     n_new_traverses = n_split_tri_list_traverses;
	     sprintf(split_label, "after %d splits",nsplits);
	     Print_Time(split_label);
	     if (tincrfiles)
	       {
		 sprintf(filename, "edges.%d",nsplits);
		 Write_Edges(filename);
		 sprintf(filename, "triangles.%d",nsplits);
		 Write_Triangles(filename);
	       }
	   }
     }
  Print_Stats();
}


//*; GET_OPTIONS  Get Input Options From Envariables 

void
Get_Options()
{
   char   *s;
   int     i;

   const numerical_opt_t *thiss;

   for (i = 0; i < n_numerical_opts; i++)
     {
	thiss = &numerical_opts[i];
	s = getenv(thiss->name);
	if (s == NULL)
	   thiss->var = thiss->defaul;
	else
	  {
	     (void) sscanf(s, "%d", &(thiss->var));
	     if (thiss->var < thiss->minval || thiss->var > thiss->maxval)
	       {
		  cout << "Error: illegal value " << thiss->var << " for" << thiss->name
		    << ", legal range=(" << thiss->minval << ',' << thiss->maxval
		    << ").\n";
		  thiss->var = thiss->defaul;
	       }
	  }
     };

   if (tinn<=1) Print_Help();
   cout << "tsplits is " << tsplits << '\n';
   Print_Time("Get_Options");
}



/* 
   INIT_VARS  Better here than in an initialization in a declaration since now I can
   restart the program w/o reloading it. */

void Init_Vars()
{
   long int tinn2 = tinn * (long int) tinn;
   short int *p;
   int     i;

   nswapq = 0;
   nswaps = 0;
   hiswapq = 0;
   max_ever_heap_size = 0;
   ndeletable = 0;
   max_allowable_num_triangles = 2 * tsplits;	// A reasonable guess.

   if (0==(points=new Point[tinn2])) 
	perror2("Allocating points array in Init_Vars.");
   if (0==(triangles=new Triangle[max_allowable_num_triangles])) 
	perror2("Allocating triangles array in Init_Vars.");

   Alloc_Dotpairs(tinn2);
   if (0==(errs = new Heap(tsplits*2)))	
	perror2("Allocating error heap in Init_Vars");
   n_elts_pushed_onto_errs_heap = 0;
   n_elts_popped_from_heap = 0;
   n_elts_deleted_from_heap = 0;
   n_elts_pushed_onto_errs_heap_swap = 0;
   n_swaps_while_popping_heap  = 0;
   n_split_tri_list_traverses = 0;
   n_new_traverses = 0;
   Print_Time("Init_Vars");

 }



// ALLOC_dotpairs  Allocate the dotted pairs array and link the freestore.

void Alloc_Dotpairs(const int n)
{
  int     i;
  dotpairs = new Dotpair[n];
  if (dotpairs==0) perror2("Unable to alloc dotpairs in alloc_dotpairs.");
  freelist.head = 0;
  for (i = 1; i < n; i++) dotpairs[i - 1].next = i;
  dotpairs[n - 1].next = eolist;
  for (i=0;i<n;i++) dotpairs[i].ptid= -9999;
}  



/** READ_DATA

1. Open stdin containing a square array of short int heights.
2. Do some validity checks on constants.
3. Read the heights into a linear array, supply Xs and Ys, and determine bounds.
*/

void
Read_Data()
{
   int     i, j;
   short int zmin, zmax;

   zmin = 32767;
   zmax = -32767;

   // Copy The Heights Data To The Points Array. 

   npoints = 0;
   for (i = 0; i < tinn; i++)
      for (j = 0; j < tinn; j++)
	{
	   points[npoints].x = i;
	   points[npoints].y = j;
	   cin >> points[npoints].z;
	if (! cin.good())
	{
	printf("At npoints=%d\n",npoints);
	perror2("Reading cin is not good.");
}

	   zmin = min(zmin, points[npoints].z);
	   zmax = max(zmax, points[npoints].z);
	   npoints++;
	}
   cout << "tinn=" << tinn << ", Z's range is [" << zmin << "," << zmax << "]\n";

   Print_Time("Read_Data");

#ifdef DEBUG
   Print_Points();
#endif
}



/* 
   CALC_ONE_EQN_AND_ERR Calculate the plane equation of triangle #t, and then find its
   maximum error, both which point and the value, and store it in the triangle.  Push 
   triangle onto errors heap.  */

void
Calc_One_Eqn_and_Err(const int t)
{
   int     i, n;
   Triangle *thist;
   COORDT     d;			// Max Error Seen So Far. 
   int     k;			// Where It Was Seen 
   Pairindex thispair;
   Point *thisp;
   COORDT     d1;
   float   eqn[3];		/* The equation of the  plane defined by the triangle's 3
  			        vertices, in the form z=x*eqn[0]+y*eqn[1]+eqn[2]. */
     
   Calc_Plane_Eqn(t, eqn);
   thist = &triangles[t];
   d = 0;
   k = -1;
   thispair = (thist -> l).head;
   while (thispair != eolist)
     {
       if (thispair<0) perror2("Negative thispair in Calc_One_Eqn_and_Err");
       thisp = &points[dotpairs[thispair].ptid];
       d1 =  fabs(eqn[0] * thisp->x + eqn[1] * thisp->y + eqn[2] - thisp->z);
       if (d1 > d)
	 {
	   d = d1;
	   k = dotpairs[thispair].ptid;
	 }
       thispair = dotpairs[thispair].next;
     }
   thist->worst_pt = k;
   thist->maxe = d;
   errs->push(t);		// Insert triangle and err into heap of errs.
}



/* 
   CALC_PLANE_EQN  Calculate equation of triangle #t from its vertices, and insert
   into the triangle.  Its form is z=ax+by+c. */

void
Calc_Plane_Eqn(const int t, float *const eqn)
{
   Triangle *thist;
   Point *a, *b, *c;	// The Vertices' Addresses (Not Their Ids). 
   int     db[3], dc[3];	// The Delta Vectors For 2 Sides. 
   float   u[3];		// Unnormalized Normal Vector. 
   float   d;			// When Eqn Is In Form Ax+By+Cz=D, This Is 'D'. 
   int     i;

   thist = &triangles[t];
   a = &points[thist->vert[0]];
   b = &points[thist->vert[1]];
   c = &points[thist->vert[2]];

   db[0] = b->x - a->x;
   db[1] = b->y - a->y;
   db[2] = b->z - a->z;

   dc[0] = c->x - a->x;
   dc[1] = c->y - a->y;
   dc[2] = c->z - a->z;
   
   Cross_Product(u, db, dc);

   if (abs(u[2]) < TOLERANCE)
     {
	if (tprint!=0) printf("Warning: Plane of triangle #%d is nearly vertical.\n", t);
	eqn[0] = 0.0;
	eqn[1] = 0.0;
	eqn[2] = 0.0;
   } else
     {
	eqn[0] = -u[0] / u[2];
	eqn[1] = -u[1] / u[2];
	d = Dot3D(a, u);
	eqn[2] = d / u[2];
     }
}




// CROSS_PRODUCT  A = B X C 

void
Cross_Product(float *const a, const float *const b, const float *const c)
{
   a[0] = b[1] * c[2] - b[2] * c[1];
   a[1] = b[2] * c[0] - b[0] * c[2];
   a[2] = b[0] * c[1] - b[1] * c[0];
}


void
Cross_Product(int *const a, const int *const b, const int *const c)
{
   a[0] = b[1] * c[2] - b[2] * c[1];
   a[1] = b[2] * c[0] - b[0] * c[2];
   a[2] = b[0] * c[1] - b[1] * c[0];
}


void
Cross_Product(float *const a, const int *const b, const int *const c)
{
   a[0] = (float) b[1] * (float) c[2] - (float) b[2] * (float) c[1];
   a[1] = (float) b[2] * (float) c[0] - (float) b[0] * (float) c[2];
   a[2] = (float) b[0] * (float) c[1] - (float) b[1] * (float) c[0];
}


// PRINT_TRIANGLE 

void
Print_Triangle(const int t)
{
   Triangle *tt;
   int     i;

   if (tprint==0) return;
   
   tt = &triangles[t];
   cout << "Triangle #" << t << ": vert: " << tt->vert[0] << ' ' << tt->vert[1] << ' ' << tt->vert[2]
     << ", adj: " <<   tt->adj[0] << ' ' << tt->adj[1] << ' ' << tt->adj[2]
     << ", max_err_over_whole_TIN=" << tt->maxe << " at " << tt->worst_pt << ", np=" << tt->l.length() << "\n";
#ifdef DEBUG
   printf("   p=");
   for (i = 0; i < tt->np; i++)
     {
	printf(" %d", tt->p[i]);
	if ((i && 15) == 0 && i > 0)
	   printf("\n   ");
     }
   if ((i == 0) || ((i && 15) != 0))
      printf("\n");
   printf("worst_pt=%d, maxe=%.3f, ", tt->worst_pt, tt->maxe);
#endif
}


// PRINT_POINTS All Of Them 

void
Print_Points()
{
   int     i;

   printf("Npoints=%d\n", npoints);

   for (i = 0; i < npoints; i++)
     {
       cout << i << ": " << points[i].x << ", " << points[i].y << ", " << points[i].z;
       if ((i & 1) == 1) cout << '\n';
     }
   cout << endl;
}



// FLAG_NULL_TRIANGLE_FOR_DELETION: Test whether triangle has no
// interior points, has zero area, and is also on the outside of the
// square.  If so, flag for later deletion. It's too messy to delete
// them now; that would interfere with putting them into the swap
// queue.  The deletable triangles are stored in an array.

// A null triangle that is not on the outside should be removed soon by
// having a diagonal swapped.

void
Flag_Null_Triangle_For_Deletion(int t)
{
   int     i, j, k;
   Triangle &tri = triangles[t];

   if (tri.l.head != eolist)
      return;
   if (tri.adj[0] >= 0 && tri.adj[1] >= 0 && tri.adj[2] >= 0)
      return;
   if (abs(Area(tri)) > TOLERANCE) return;

   if (tprint!=0) printf("Flagging triangle #%d for deletion since it has no points, has zero area, and is on the outside.\n", t);

   if (ndeletable >= MDELETABLE)
     {
	printf("ERROR: ndeletable too large.\n");
	exit(1);
     }
   deletable_triangle[ndeletable++] = t;
}



// DELETE_FLAGGED_TRIANGLES 

void
Delete_Flagged_Triangles()
{
  int     i, j, t;
  Triangle *thist;
  
  for (i = 0; i < ndeletable; i++)
    {
      t = deletable_triangle[i];
      thist = &triangles[t];
      if ((thist->l).head != eolist)
	cout << "Warning, in Delete_Flagged_Triangles, deleting triangle #" << t
	     << ", but it has a nonnull point list.\n";
      
      /* Any triangle that was adjacent to this triangle is now adjacent to the
	 outside. */
      
      Update_Adjacency(t, -1);
      
      if (t < ntriangles - 1)
	{			// Copy Last Triangle Into This Position. 
	    Copy_Triangle(t, ntriangles -1);
	  
	  /* Any triangle that was adjacent to #ntriangles-1 is now adjacent to
	     #t. */
	    Update_Adjacency(ntriangles - 1, t);
	}
      ntriangles--;
      if (tprint>0)
	cout << "Deleting triangle #" << t << ", nt now=" << ntriangles << '\n';
    }
  ndeletable = 0;
}



/* 
   UPDATE_ADJACENCY  The triangle in triangles[oldt] has been moved to
   triangles[newt].  Now update all references to it.  For every triangle that was
   adjacent to #oldt, make it adjacent to #newt in that entry.  Also change every
   instance of #oldt to #newt in swapq.  If newt= -1 don't update swapq since there
   should be no oldt's in swapq. */

void
Update_Adjacency(const int oldt, const int newt)
{
   int     i, j, t;

   for (i = 0; i < 3; i++)
     {
	t = triangles[oldt].adj[i];
	if (t >= 0)
	  {			// Triangle T Is Adjacent.  Update It. 
	     for (j = 0; j < 3; j++)
		if (triangles[t].adj[j] == oldt)
		  {
		     triangles[t].adj[j] = newt;
#ifdef DEBUG
		     printf("Updating adjacency of triangle %d from %d to %d.\n", t, oldt, newt);
#endif
		     goto foundadj;
		  }
	     printf("ERROR: In Update_Adjacency, triangle #%d is adjacent to triangle #%d, "
		    "but the reverse is not true.\n", oldt, t);
	   foundadj:;
	  }
     }

   if (newt >= 0)
      for (i = 0; i < nswapq; i++)
	 for (j = 0; j < 2; j++)
	    if (swapq[i][j] == oldt)
	      {
		 swapq[i][j] = newt;
		 if (tprint!=0) printf("In update_adjacencies, swapq[%d][%d]=oldt=%d changed to %d.\n",
			i, j, oldt, newt);
	      }
}


/* 
   ADD_TO_SWAPQ Add a triangle pair, defining a new edge, and a new diagonal of a
   quadrilateral, to the swap queue, to be considered for swapping.  The inputs are
   one triangle number and the number of its adj entry defining the other triangle.
   Don't add an edge if it's already in the queue. */

void
Add_to_Swapq(const int t, const int a)
{
   int     i;

   if (tswap==0) return;		// Maybe user doesn't want any swaps.
   
   i = triangles[t].adj[a];
   if (i < 0) return;			// Neighbor Is The Outside. 
   if (nswapq >= MSWAPQ)
     {
	printf("ERROR: MSWAPQ=%d is too small.\n", MSWAPQ);
	exit(1);
     }
   swapq[nswapq][0] = min(t, i);
   swapq[nswapq][1] = max(t, i);

   // Is This Edge Already There? 

   for (i = 0; i < nswapq; i++)
      if (swapq[i][0] == swapq[nswapq][0] && swapq[i][1] == swapq[nswapq][1])
	{
#ifdef DEBUG
	   printf("Not adding (%d,%d) to swapq since it duplicates entry #%d\n",
		  swapq[nswapq][0], swapq[nswapq][1], i);
#endif
	   goto dupswap;
	}
#ifdef DEBUG
   printf("Adding (%d,%d) to swapq[%d].\n", swapq[nswapq][0], swapq[nswapq][1], nswapq);
#endif
   nswapq++;
   hiswapq = max(hiswapq, nswapq);
 dupswap:;
}


/* 
   WRITE_PGM Write the data about which triangle each point is in as a portable
   greymap (PGM) file so that it can be displayed for debugging purposes, by, e.g.,
   xv, which can pseudocolor the greyscales.  Write_Pgm works only when the points
   are in an array.  Assign the triangle numbers as their colors. Make the vertices
   white, coded by ntriangles in c. */

void
Write_Pgm()
{
   int     i, j, k, itr, x, y;
   int     c[tinn][tinn];	// will this work? 
   FILE   *f;

   /* Zero c, which shouldn't be necessary since every entry should be set later,
      but be careful. */

   for (i = 0; i < tinn; i++)
      for (j = 0; j < tinn; j++)
	 c[i][j] = 0;

   for (itr = 0; itr < ntriangles; itr++)
     {
	for (j = 0; j < 3; j++)
	  {			// Color The Vertices. 
	     k = triangles[itr].vert[j];
	     x = points[k].x;
	     y = points[k].y;
	     if (x < 0 || x >= tinn || y < 0 || y >= tinn)
		printf("ERROR: illegal (x,y)=(%d,%d) in Write_Pgm\n", x, y);
	     else
		c[x][y] = 0;
	  }
	Pairindex thispair;
	thispair = triangles[itr].l.head;
	while (thispair != eolist)
	  {			// Color The Points In The Triangles. 
	    if (thispair<0) perror2("Negative thispair in Write_Pgm.");
	    k = dotpairs[thispair].ptid;
	    thispair=dotpairs[thispair].next;
	    x = points[k].x;
	    y = points[k].y;
	    if (x < 0 || x >= tinn || y < 0 || y >= tinn)
	      printf("ERROR: illegal (x,y)=(%d,%d) in Write_Pgm\n", x, y);
	    else
	      c[x][y] = itr;
	  }
     }

   ofstream tout("t.pgm");
   if (!tout) perror2("Can't open t.pgm file for writing.");
   tout << "P2\n#t.pgm, written by tin.cc\n" << tinn << ' ' << tinn 
      << endl << ntriangles << endl;

   k = 0;			// Counter For When To Insert A Newline. 
   for (i = 0; i < tinn; i++)
      for (j = 0; j < tinn; j++)
	{
	  tout << setw(8) << c[i][j];
	  if (((++k) & 7) == 0) tout << endl;
	}
   if ((k & 7) != 0) tout << endl;
}


/* WRITE_EDGES Write the triangles' edges to a file, with each edge:

  x1  y1  z1
  x2  y2  z2
  <blank line>

 Altho each edge (except on the outside) occurs in 2 triangles, write it
 only once.  Do this by writing an edge iff the adjacent triangle has a
 smaller number, which includes -1 for the outside. */

void Write_Edges(const char *const filename)
{
  int itr;
  int j;
  int a,b;
  ofstream to(filename);

  if (!to) perror2("Can't open edges file for writing.");

  to << "# Edges in the triangulation:" << endl
     << "# x0 y0 z0" << endl 
     << "# x1 y1 z1" << endl
     << "# (blank line)" << endl
     << "# ..." << endl ;
  
  for (itr=0; itr<ntriangles; itr++)
    for (j=0; j<3; j++)
      {
	a = triangles[itr].vert[j];
	b = triangles[itr].vert[plus1[j]];
	if (itr > triangles[itr].adj[minus1[j]]) 
	  to << setw(5) << points[a].x << setw(5) << points[a].y << setw(5) << points[a].z << endl
	     << setw(5) << points[b].x << setw(5) << points[b].y << setw(5) << points[b].z << endl << endl;
      }	  
}



/* WRITE_TRIANGLES Write the triangles to a file, with each triangle:

  x1  y1  z1
  x2  y2  z2
  x3  y3  z3
  x1  y1  z1
  <blank line>
  */

void Write_Triangles(const char *const filename)
{
  int itr;
  int j;
  int a,b;
  ofstream to(filename);

  if (!to) perror2("Can't open triangles file for writing.");

  to << "# Triangles in the triangulation:" << endl
     << "# x0 y0 z0" <<endl
     << "# x1 y1 z1" <<endl
     << "# x2 y2 z2" <<endl
     << "# x0 y0 z0" <<endl
     << "# (blank line)" <<endl
     << "# ..." << endl;

  for (itr=0; itr<ntriangles; itr++)
    {
      for (j=0; j<3; j++)
	{
	  a = triangles[itr].vert[j];
	  to << setw(5) << points[a].x << setw(5) << points[a].y 
	     << setw(5) << points[a].z << endl;
	}	  
      a = triangles[itr].vert[0];
      to << setw(5) << points[a].x << setw(5) << points[a].y 
	 << setw(5) << points[a].z << endl << endl;
    }
}




// WRITE_TRIANGLE_CENTERS Write the triangles' centers to a file in
// the form of gnuplot set label commands.  For triangle #7:

// set label 8 "7"   x  y  z


void Write_Triangle_Centers(const char *const filename)
{
  int itr;
  int b;
  float a;
  ofstream to(filename);

  if (!to) perror2("Can't open triangle centers file for writing.");

  to << "# tin triangle centers.   do this in gnuplot:" << endl
     << "# load \"centers.all\"" << endl;

  for (itr=0; itr<ntriangles; itr++)
      {
	to << "set label " << setw(7) << itr+1 << ' ' << '"' << setw(7) << itr << '"'; 
	  a = (points[triangles[itr].vert[0]].x + points[triangles[itr].vert[1]].x +
	       points[triangles[itr].vert[2]].x) * 0.33333;
	  to << " at " << setw(7) << a;
	  a = (points[triangles[itr].vert[0]].y + points[triangles[itr].vert[1]].y + 
	       points[triangles[itr].vert[2]].y) * 0.33333;
	  to << ',' << setw(7) << a;
	  a = (points[triangles[itr].vert[0]].z + points[triangles[itr].vert[1]].z +
	       points[triangles[itr].vert[2]].z) * 0.33333;
	  to << ',' << setw(7) << a << " center\n";
    }
}




void Print_Stats()
{
  Print_Time2("cumulative processing (not i/o)");
  cout << "[Tin finished.]\n"
       << "# rows of input points (envar TINN) =                " << tinn << '\n'
       << "# input points =                                     " << npoints << '\n'
       << "Max # splits allocated (envar TSPLITS) input as =    " << tsplits << '\n'
       << "# splits =                                           " << nsplits << '\n'
       << "# output triangles =                                 " << ntriangles << '\n'
       << "When are quadrilateral diagonals swapped?            " ;

  switch(tswap)
  {   
  case 0: cout << "never"; break;
  case 1: cout << "Delauney"; break;
  case 2: cout << "min weight"; break;
  case 3: cout << "min max error"; break;
  default: cout<< "illegal value of tswap";
  }
  cout << endl;

  cout << "# times quadrilateral diagonal swapped =             " << nswaps << '\n'
       << "Hi-water size for edge swap queue =                  " <<  hiswapq << '\n'
       << "# times elements pushed onto errs heap =             " <<  n_elts_pushed_onto_errs_heap << '\n'
       << "# times elements popped from errs =                  " << n_elts_popped_from_heap << '\n'
       << "# times elements deleted from errs =                 " << n_elts_deleted_from_heap << '\n'
       << "# times elements swapped when pushing =              " << n_elts_pushed_onto_errs_heap_swap << '\n'
       << "# times elements swapped when popping =              " << n_swaps_while_popping_heap << '\n'
       << "Hi-water size for heap =                             " << max_ever_heap_size << '\n'
       ;
  Write_Edges("edges.all");
  Write_Triangles("triangles.all");
  Write_Triangle_Centers("centers.all");
  Write_Pgm();
  Print_Time("at end");
}


// HEAP_ELT_COMPARE Call by qsort called by heap push.  This is defined so as to
// sort in DECREASING order.  If two heap elements have the same error, then
// don't worry about order.

int Heap_Elt_Compare(const void *a, const void *b)
{
  int d= (int) ( triangles[((const Heap_elt *)b)->tri].maxe - 
		 triangles[((const Heap_elt *)a)->tri].maxe);
  return d;
}




/** 
   MAKE_FIRST_2_TRIANGLES  This assumes that the input heights are in a square array
   with y varying faster than x. 


          ...........
          .2      1..
          .       .2.
          .  T1  .  .
          .     .   .
          .    .    .
          .   .     .
          .  .  T0  .
          .0.       .
          ..0      1.
          ...........
*/


void
Make_First_2_Triangles()
{
   Triangle *t0, *t1;		/* Pointers to the triangles, which are really in
				   array 'triangles'. */
   int     i, t;


   if (max_allowable_num_triangles < 2)
     {
	printf("max_allowable_num_triangles=%d is ridiculously small.\n", max_allowable_num_triangles);
	exit(1);
     }
   t0 = &triangles[0];
   t1 = &triangles[1];

   t0->vert[0] = 0;
   t0->vert[1] = tinn * (tinn - 1);
   t0->vert[2] = tinn * tinn - 1;
   t1->vert[0] = 0;
   t1->vert[1] = tinn * tinn - 1;
   t1->vert[2] = tinn - 1;
   t0->adj[0] = -1;
   t0->adj[1] = 1;
   t0->adj[2] = -1;
   t1->adj[0] = -1;
   t1->adj[1] = -1;
   t1->adj[2] = 0;

   /* 
      Allocate the points to the 2 triangles based on their indices since that
      determines their X and Y coordinates.  The vertices of a triangle are not in
      any triangle. */

   for (i = 0; i < 2; i++)
     triangles[i].l = eolist;

   for (i = 0; i < npoints; i++)
     if (i == 0 || i == tinn -1 || i == tinn * tinn - 1 || i == tinn * (tinn - 1)); 
      else if (points[i].x > points[i].y)
	{
	  triangles[0].l.push(i);
      } else
	{
	  triangles[1].l.push(i);
	}

   ntriangles = 2;
   Calc_One_Eqn_and_Err(0);
   Calc_One_Eqn_and_Err(1);
   Print_Triangle(0);
   Print_Triangle(1);

   Print_Time("make_1st_2_triangles");
}



/**

SPLIT_TRIANGLE Split the worst triangle, calculate the new
triangles' equations etc, and update max_err_over_whole_TIN etc.
How the 3 new triangles' vertices correspond to the old triangle,
and the signs of the new edge equations. Return 0 if this was
successful, or 1 if no triangle could be split.  This might be
because no triangle has any error.



                                2
                                .
                               ...
                              .2.2.
                             .  .  .
                            .   .   .
                           .    .    .
                        1 .    -.+    . 0
                         .      .      .
new triangle #ntriangles.      1.0      .       new triangle #worst_triangle
                       .      ....       .
                      .    +..  2 ..-     .
                     .    ..-      +..     .
                    .   ..            ..    .
                   . 0..                ..   .
                  . ..                    ..1 .
                 ... 0                     1.. .
                .................................
               0                2                1
                     new triangle #ntriangles+1 */


int Split_Triangle()
{
   int     i, j, k;
   int     newt;		// new triangle 
   Triangle *thist;
   float  line[3][3];		/* 2D equations of the 3 new edges: line[i][0]*x +
				   line[i][1]*y + line[i][2] = 0 */
   int     newvert;		/* ID of worst point of worst triangle; the new
				   vertex.  */
   Point *thispt;		// Address Of Coords Of New Point. 
   int     nnewp[3];		// Number Of Points In Each New Triangle. 
   int     onedge[3];		// 1=> the new point is on edge #i of the triangle. 
   float  d;			// Distance Of Current Point From Line. 

   if (tprint>0) cout << "Trying to split the worst triangle.\n";

   if (worst_triangle < 0 || max_err_over_whole_TIN == 0.0 || triangles[worst_triangle].worst_pt < 0)
     {
       cout << "No triangle splittable.\n";
       return 1;
     }
   thist = &triangles[worst_triangle];

   // Allocate 2 New Triangles And Split The Worst Triangle. 

   if (ntriangles + 2 > max_allowable_num_triangles)
     {
	printf("ERROR: max_allowable_num_triangles=%d is too small.\n", max_allowable_num_triangles);
	exit(1);
     }
   newvert = thist->worst_pt;
   if (tprint!=0) cout << "Splitting triangle " << worst_triangle << ", with max_err_over_whole_TIN "
		       << max_err_over_whole_TIN << ", at point " << newvert << '\n';

   thispt = &points[newvert];

   /* 
      Find 2D equations of the 3 new edges.  Must do before updating the old
      triangle.  This doesn't overflow because, altho points is short int, line is
      float. */

   for (i = 0; i < 3; i++)
     {
	line[i][0] = points[thist->vert[i]].y - thispt->y;
	line[i][1] = -points[thist->vert[i]].x + thispt->x;
	line[i][2] = -thispt->x * line[i][0] - thispt->y * line[i][1];
     }

#ifdef DEBUG
   printf("Line equations of new edges:\n");
   for (i = 0; i < 3; i++)
      printf("%d: %.3f*x%+.3f*y%+.3f=0.\n", i, line[i][0], line[i][1], line[i][2]);
#endif

   // Test Which Side Of Line L Point P Is On. 

#define TESTLP(l,p) (line[l][0]*points[p].x+line[l][1]*points[p].y+line[l][2])


   // Update The Adjacencies Of The Triangles Adjacent To Worst_Triangle. 

   for (i = 0; i < 3; i++)
     {
	k = triangles[worst_triangle].adj[i];
	if (k >= 0)
	  {
	     // Triangle K Is Adjacent.  Update It. 
	     switch (i)
	       {
	       case 0:
		  newt = worst_triangle;
		  break;
	       case 1:
		  newt = ntriangles;
		  break;
	       case 2:
		  newt = ntriangles + 1;
		  break;
	       default:
		  perror2("Shouldn't have hit default in switch in Split_Triangle.");
		  exit(1);
	       }
	     for (j = 0; j < 3; j++)
		if (triangles[k].adj[j] == worst_triangle)
		  {
		     triangles[k].adj[j] = newt;
#ifdef DEBUG
		     printf("Updating adjacency of triangle %d from %d to %d.\n", k, worst_triangle, newt);
#endif
		     goto foundadj2;
		  }
	     printf("ERROR: In Split_Triangle, triangle #%d is adjacent to triangle #%d, but the reverse is not true.\n",
		    worst_triangle, k);
	   foundadj2:;
	  }
     }

   // Create The 3 New Triangles.  One Is The Old Triangle, Updated. 

   triangles[ntriangles].vert[0] = thist->vert[0];
   triangles[ntriangles].vert[1] = newvert;
   triangles[ntriangles].vert[2] = thist->vert[2];

   triangles[ntriangles + 1].vert[0] = thist->vert[0];
   triangles[ntriangles + 1].vert[1] = thist->vert[1];
   triangles[ntriangles + 1].vert[2] = newvert;

   thist->vert[0] = newvert;	// Do This After The Previous Two Triangles. 
   // thist->vert[1], [2] are the same. 

   triangles[ntriangles].adj[0] = worst_triangle;
   triangles[ntriangles].adj[1] = thist->adj[1];
   triangles[ntriangles].adj[2] = ntriangles + 1;

   triangles[ntriangles + 1].adj[0] = worst_triangle;
   triangles[ntriangles + 1].adj[1] = ntriangles;
   triangles[ntriangles + 1].adj[2] = thist->adj[2];

   // thist->adj[0] is the same. 
   thist->adj[1] = ntriangles;
   thist->adj[2] = ntriangles + 1;

   if (fabs(TESTLP(0, thist->vert[1])) < TOLERANCE)
      onedge[2] = 1;
   else
      onedge[2] = 0;
   if (fabs(TESTLP(0, thist->vert[2])) < TOLERANCE)
      onedge[1] = 1;
   else
      onedge[1] = 0;
   if (fabs(TESTLP(1, thist->vert[2])) < TOLERANCE)
      onedge[0] = 1;
   else
      onedge[0] = 0;

   for (i = 0; i < 3; i++)
      if (onedge[i] > 0)
	 if (tprint!=0) printf("New vertex is on edge #%d of triangle.\n", i);

   /* 
      Allocate each point from the old triangle, except the new vertex, to one of
      the three new triangles.  Points that are on a new edge can be put into either 
      adjacent triangle unless one adjacent triangle is degenerate.  In that case,
      put the point in the other triangle.  The problem that this avoids is as 
      follows.  A degenerate triangle is probably a vertical plane; I can't
      calculate a plane equation, so I can't calculate errors of points from it
      (or rather they are infinite). */

   // In terms of the former version: 0: Worst, 1: Ntriangles, 2: Ntriangles+1 

   Linked_list list0, list1, list2;
   list0 = list1 = list2 = eolist;

   while (thist->l.head != eolist) 
     {
       n_split_tri_list_traverses++;
       j = thist->l.pop();
       if (j < 0) perror2("Illegal point popped in Split_Triangle.");
       if (j == newvert) ;
       else
	 {
	   d = TESTLP(0, j);
	   if (d > TOLERANCE || (d >= 0 && onedge[1] == 0))
	     {
	       d = TESTLP(2, j);
	       if (d > TOLERANCE || (d >= 0 && onedge[0] == 0)) list0.push(j);
	       else list1.push(j);
	     } else
	       {
		 d = TESTLP(1, j);
		 if (d > TOLERANCE || (d >= 0 && onedge[2] == 0)) list2.push(j);
		 else list0.push(j);
	       }
	 }
     }

   triangles[worst_triangle].l = list0;
   triangles[ntriangles].l = list1;
   triangles[ntriangles+1].l = list2;
   
//already deleted when popped:   errs->del(triangles[worst_triangle].heap_loc);

   Calc_One_Eqn_and_Err(worst_triangle);
   Calc_One_Eqn_and_Err(ntriangles);
   Calc_One_Eqn_and_Err(ntriangles + 1);

   /* 
      Add edges between the three new triangles and the adjacent old triangles to
      the list of quadrilateral diagonals to be tested to swapping. */

   Add_to_Swapq(worst_triangle, 0);
   Add_to_Swapq(ntriangles, 1);
   Add_to_Swapq(ntriangles + 1, 2);

   ntriangles += 2;

   Print_Triangle(worst_triangle);
   Print_Triangle(ntriangles - 2);
   Print_Triangle(ntriangles - 1);

/* Check triangles last to first, else if there are 2 null triangles, then deleting
   the first will invalidate the deletable_triangle array. */

   Flag_Null_Triangle_For_Deletion(ntriangles - 1);
   Flag_Null_Triangle_For_Deletion(ntriangles - 2);
   Flag_Null_Triangle_For_Deletion(worst_triangle);

   Delete_Flagged_Triangles();

   return 0;

}




/**
  PROC_SWAPQ  Test diagonals in swapq and maybe swap them.
*/


void
Proc_Swapq()
{
   int     iswap, i, j, k, ixy, itr, n, t, ta;
   int     v0, v1;		// ids of 2 current vertices defining an edge. 
   int     t0, t1;		// the ids of the 2 adjacent triangles. 
   Triangle *(thist[2]);	/* pointers to the 2 triangles.  'thist' is an array
			       of 2 pointers to 'triangle'. */
   int     a[2];		/* which adj entry in each triangle points to the
				other. */
   int doswap;			// Do we swap this diagonal?  1=yes 0=no.
   float   oldcos[2][3], newcos[2][3];	/* cosines of the angles of the existing,
				       and proposed new, triangles. */
   float   oldmax, newmax;	// Maxima Of The Old And New Cosines. 
   COORDT  oldedge[2][3][2];	/* The edge vectors of the two current triangles:
			       oldedge[which triangle][which edge of it][x or y]. */
   COORDT  newedge[2][3][2];	// Edge Vectors Of Proposed New Triangles. 
   int     newv[2][3], newa[2][3];	// The New Vert And Adj. 
   COORDT  line[3];		/* Equation of new diagonal:
			       line[0]*x+line[1]*y+line[2]=0. */
   Point *thisp;		// Pointer To Current Point, In Array 'Points'. 
   int    *(oldp[2]);		/* Copy of pointers to p arrays for old triangles,
				while I'm creating new p arrays. */
   COORDT det;			// Determinant of a pair  of edges in one of the
				// new triangles.
   COORDT mindet;		// Smallest determinant.
   float  oldsqlen, newsqlen;   // Squared lengths of the old and new diagonals.

   oldmax = newmax = -99999e9; // stop a compiler warning
   /* 
      Repeat for each candidate diagonal.  WARNING:  The body of this loop may add
      diagonals to the q, and so increase nswapq. */

   for (iswap = 0; iswap < nswapq; iswap++)
     {
        doswap = 0;
	t0 = swapq[iswap][0];
	t1 = swapq[iswap][1];
	thist[0] = &triangles[t0];
	thist[1] = &triangles[t1];

	// For Each Of The 2 Triangles, Find Which Adj Points To The Other. 

	for (itr = 0; itr < 2; itr++)
	  {			// For Each Of The Diagonal'S Adjacent Triangles 
	     for (j = 0; j < 3; j++)	/* Check all of this triangle's adjacent
					   triangles to find the match. */
		if (thist[itr]->adj[j] == swapq[iswap][1 - itr])
		  {
		     a[itr] = j;
		     goto foundswap;
		  }
	     printf("ERROR: in Proc_Swapq, cannot find a[%d] for swapq[%d][%d]=%d\n", itr, iswap, 1 - itr);
	     exit(1);
	   foundswap:;
	  }

	/* 
	   Find the edge vectors of each old triangle.  Make each edge be opposite
	   the corresponding vertex. */

	for (itr = 0; itr < 2; itr++)
	   for (j = 0; j < 3; j++)
	     {
		v1 = thist[itr]->vert[plus1[j]];
		v0 = thist[itr]->vert[minus1[j]];
		oldedge[itr][j][0] = points[v1].x - points[v0].x;
		oldedge[itr][j][1] = points[v1].y - points[v0].y;
	     }

	/* 
	   Calculate the new vert into a new array, then later copy them back if we
	   do the swap. */

	newv[0][0] = thist[0]->vert[minus1[a[0]]];
	newv[0][1] = thist[0]->vert[a[0]];
	newv[0][2] = thist[1]->vert[a[1]];
	newv[1][0] = thist[0]->vert[a[0]];
	newv[1][1] = thist[0]->vert[plus1[a[0]]];
	newv[1][2] = thist[1]->vert[a[1]];

	/* 
	   Calculate edge vectors of the hypothetical new triangles. */

	for (itr = 0; itr < 2; itr++)
	   for (j = 0; j < 3; j++)
	     {
		v1 = newv[itr][plus1[j]];
		v0 = newv[itr][minus1[j]];
		newedge[itr][j][0] = points[v1].x - points[v0].x;
		newedge[itr][j][1] = points[v1].y - points[v0].y;
	     }

	// Now, handle the different swap criteria.

	switch (tswap)  {
	case 0: break;          // never swap
	case 1: 		// Delauney
	  /* 
	     Test, and swap diagonals if necessary to maintain a Delauney
	     triangulation, i.e., to maximize the minimum angle in the 2 triangles. To 
	     save a trig evaluation, I'll minimize the maximum cosine of the angles.
	     
	     Find the angle cosines.  The minus is needed because each edge vector is
	     coming from a different vertex. */

	  oldmax = -100.0;
	  for (itr = 0; itr < 2; itr++)
	    for (j = 0; j < 3; j++)
	     {
		oldcos[itr][j] =
		   -Dot2D(oldedge[itr][minus1[j]], oldedge[itr][plus1[j]]) /
		   sqrt(Dot2D(oldedge[itr][minus1[j]], oldedge[itr][minus1[j]]) *
			Dot2D(oldedge[itr][plus1[j]], oldedge[itr][plus1[j]]));
		if (oldcos[itr][j] > oldmax)
		   oldmax = oldcos[itr][j];
	     }

	newmax = -100.;
	mindet = 9999;	// Precise value unimportant; need to know only if negative.

	for (itr = 0; itr < 2; itr++)
	    {
		for (j = 0; j < 3; j++)
		    {			// Find The Angle Cosines 
			newcos[itr][j] =
			    -Dot2D(newedge[itr][j], newedge[itr][plus1[j]]) /
			    sqrt(Dot2D(newedge[itr][j], newedge[itr][j]) *
				 Dot2D(newedge[itr][plus1[j]], newedge[itr][plus1[j]]));
			if (newcos[itr][j] > newmax)
			    newmax = newcos[itr][j];
		    }
// Find the determinants of 3 pairs of edges for each of the 2 new
// triangles, in order to check whether a new triangle is oriented
// backwards.  This will happen if the proposed new diagonal lies outside
// the triangles.  (You might want to check the adjacent new angles are
// acute, but this is too strict.)  It's sufficient to do the 2 pairs
// representing the new angles, but determining which two is a little more
// complicated.

			det = Det2D(newedge[itr][0], newedge[itr][1]);
			Mineq(mindet, det);
			det = Det2D(newedge[itr][1], newedge[itr][2]);
			Mineq(mindet, det);
			det = Det2D(newedge[itr][2], newedge[itr][0]);
			Mineq(mindet, det);
	    }

#ifdef DEBUG
	printf("oldcos= %.3f %.3f %.3f %.3f %.3f %.3f\nnewcos=%.3f %.3f %.3f %.3f %.3f %.3f\n",
	       oldcos[0][0], oldcos[0][1], oldcos[0][2], oldcos[1][0], oldcos[1][1], oldcos[1][2],
	       newcos[0][0], newcos[0][1], newcos[0][2], newcos[1][0], newcos[1][1], newcos[1][2]);
#endif

// If the min new angle is larger than the min old angle, and the new
// triangles have positive angles, then swap.

if (newmax < oldmax && mindet > 0) doswap=1;  
 break;

	case 2:   	// min weight

	  oldsqlen = Square(oldedge[0][0][0])+Square(oldedge[0][0][1]);
	  newsqlen = Square(newedge[0][0][0])+Square(newedge[0][0][1]);
	  if (newsqlen < oldsqlen) doswap=1;
	  break;

	case 3:         // min max error
	  // not implemented yet
	default: cout << "Illegal tswap in proc_swapq " << tswap <<endl; exit(1);
	}


        if (doswap)
	  {
	     if (tprint!=0) 
                 cout << "Swap #" << nswaps << ": Swapping triangle #" << t0
                      << " (heap_loc=" << triangles[t0].heap_loc 
                      << ") and triangle #" <<  t1 << " (heap_loc=" 
                      << triangles[t1].heap_loc << ")"<< endl;
	     nswaps++;


	    /**
             v[a[0]] means thist[0]->vert[a[0]] etc

                                OLD                                   NEW

                                        v[minus1[a[1]]]            2
                    v[a[1]]......................                ......................1
                           .                   ..v[plus1[a[0]]]  ..                   .
                           .                  . .               2. .                  .
                           .                 .  .                .  .                 .
                           .                .   .                .   .                .
                           .               .    .                .    .               .
                           .      T1      .     .                .     .       T1     .
                           .             .      .                .      .             .
                           .            .       .                .       .            .
                           .           .        .                .        .           .
                           .          .         .                .         .          .
                           .         .          .                .          .         .
                           .        .           .                .           .        .
                           .       .            .                .            .       .
                           .      .             .                .             .      .
                           .     .     T0       .                .       T0     .     .
                           .    .               .                .               .    .
                           .   .                .                .                .   .
                           .  .                 .                .                 .  .
                           . .                  .                .                  . .
                           ..                   .                .                   ..
             v[plus1[a[1]]]......................                ...................... 0
                           v[minus1[a[0]]]     v[a[0]]          0                   1


            */

	     /* The following errors mean that my internal data structure is
	        invalid. */

	     if (thist[0]->vert[plus1[a[0]]] != thist[1]->vert[minus1[a[1]]])
	       {
		  printf("ERROR: in Proc_Swapq, thist[0]->vert[plus1[a[0]]]!=thist[1]->vert[minus1[a[1]]]\n");
		  exit(1);
	       }
	     if (thist[0]->vert[minus1[a[0]]] != thist[1]->vert[plus1[a[1]]])
	       {
		  printf("ERROR: in Proc_Swapq, thist[0]->vert[minus1[a[0]]]!=thist[1]->vert[plus1[a[1]]]\n");
		  exit(1);
	       }
	     /* 
	        Calculate the new adj into a new array, then copy vert and adj back. 
	      */

	     newa[0][0] = t1;
	     newa[0][1] = thist[1]->adj[minus1[a[1]]];
	     newa[0][2] = thist[0]->adj[plus1[a[0]]];
	     newa[1][0] = thist[1]->adj[plus1[a[1]]];
	     newa[1][1] = t0;
	     newa[1][2] = thist[0]->adj[minus1[a[0]]];

	     /* 
	        Flipping the diagonal may change swapq entries referring to these
	        triangles.  Update them. */

	     for (i = iswap + 1; i < nswapq; i++)
		if (swapq[i][0] == t1 && swapq[i][1] == thist[1]->adj[minus1[a[1]]])
		  {
		     swapq[i][0] = t0;
		     if (tprint!=0) 
		       printf("Updating swapq[%d][%d] from %d to %d because of diagonal swap.\n", i, 0, t1, t0);
		} else if (swapq[i][1] == t1 && swapq[i][0] == thist[1]->adj[minus1[a[1]]])
		  {
		     swapq[i][1] = t0;
		     if (tprint!=0) 
		       printf("Updating swapq[%d][%d] from %d to %d because of diagonal swap.\n", i, 1, t1, t0);
		} else if (swapq[i][0] == t0 && swapq[i][1] == thist[0]->adj[minus1[a[0]]])
		  {
		     swapq[i][0] = t1;
		     if (tprint!=0) 
		       printf("Updating swapq[%d][%d] from %d to %d because of diagonal swap.\n", i, 0, t0, t1);
		} else if (swapq[i][1] == t0 && swapq[i][0] == thist[0]->adj[minus1[a[0]]])
		  {
		     swapq[i][1] = t1;
		     if (tprint!=0) 
		       printf("Updating swapq[%d][%d] from %d to %d because of diagonal swap.\n", i, 1, t0, t1);
		  }
	     // Now Update The Vert And Adj Entries. 

	     for (itr = 0; itr < 2; itr++)
		for (j = 0; j < 3; j++)
		  {
		     thist[itr]->vert[j] = newv[itr][j];
		     thist[itr]->adj[j] = newa[itr][j];
		  }

	     /* 
	        Update the adjacencies of the two new triangles with their
	        neighbors. Don't  update their adjacencies with each other; that's
	        the j==itr test.  It uses the fact that adj[0] of t0 is t1 

	        etc. */

	     for (itr = 0; itr < 2; itr++)
	       {
		  t = swapq[iswap][itr];
		  for (j = 0; j < 3; j++)
		    {
		       if (j == itr)
			  continue;
		       ta = thist[itr]->adj[j];
		       if (ta < 0)
			  continue;
		       for (k = 0; k < 3; k++)
			 {
			    if (triangles[ta].adj[k] == t0 || triangles[ta].adj[k] == t1)
			      {
				 triangles[ta].adj[k] = t;
				 goto found1adj;
			      }
			 }
		       printf("ERROR: When updating 1 adjacency, triangles[%d].adj[%d]=%d is not matched by the other triangle.\n", t, j, ta);
		       exit(1);
		     found1adj:;
		    }
	       }

	     /* 
	        Find the equation of the new diagonal. The positive side is triangle 
	        t1; the negative side t0. */

	     line[0] = points[thist[0]->vert[2]].y - points[thist[0]->vert[1]].y;
	     line[1] = points[thist[0]->vert[1]].x - points[thist[0]->vert[2]].x;
	     line[2] = -points[thist[0]->vert[2]].x * line[0] - points[thist[0]->vert[2]].y * line[1];
#ifdef DEBUG
	     printf("New diagonal is %.3f*x%+.3f*y%+.3f=0.\n", line[0], line[1], line[2]);
#endif

// Reassign Each Point From The Two Old Triangles To The Proper New Triangle.  

	     Linked_list linked_list0, linked_list1;
	     int j;
	     
	     linked_list0 = eolist;
	     linked_list1 = eolist;
	     
	     for (itr = 0; itr < 2; itr++)
	       {
		 while (thist[itr]->l.head != eolist) 
		   {
		      j = thist[itr]->l.pop();
		       thisp = &points[j];
		       if (line[0] * thisp->x + line[1] * thisp->y + line[2] > 0)
			 {
			   linked_list1.push(j);
		       } else
			 {
			   linked_list0.push(j);
			 }
		    }
	       }

	     triangles[t0].l = linked_list0;
	     triangles[t1].l = linked_list1;
	     
	     errs->del(triangles[t0].heap_loc);    // Delete old triangles from errs heap.
	     errs->del(triangles[t1].heap_loc);

	     Calc_One_Eqn_and_Err(t0);  
	     Calc_One_Eqn_and_Err(t1);

	     Print_Triangle(t0);
	     Print_Triangle(t1);

	     /* 
	        Add new neighboring edges of the two triangles to swapq.  Do not add 
	        the diagonal between them. */

	     for (itr = 0; itr < 2; itr++)
		for (j = 0; j < 3; j++)
		   if (thist[itr]->adj[j] != swapq[iswap][1 - itr])
		      Add_to_Swapq(swapq[iswap][itr], j);
	} else
	    if (tprint>0) 
		cout << "Not swapping triangle #" << t0 << " and triangle #" << t1
		     << ".  Oldmax= " << oldmax << ", newmax = " << newmax
		     << ", mindet= " << mindet << '\n';

     }				// End of processing this swapq entry. 
   nswapq = 0;
}




// PRINT_HELP  Copy the file tin.help to stdout and exit.

void Print_Help()
{
  char *argv[]={"bincat","tin.help",NULL};
  execve("/bin/cat",argv,NULL);
}




// HEAP:PUSH
// Push given triangle onto the heap, if its error > 0.

void Heap::push(const int t)
{
 int father, current;

 if (triangles[t].maxe == 0.0)
 {  if (tprint>0) cout << "Not pushing triangle #" << t << " since its err=0." << endl;
    triangles[t].heap_loc = -1;
    return; }

 if (heap_size >= alloc_size) perror2("Heap overflow.");
 current = ++heap_size;
 max_ever_heap_size = max(max_ever_heap_size, heap_size);
 elts[heap_size].tri = t;
 triangles[t].heap_loc = heap_size;
 n_elts_pushed_onto_errs_heap++;
 if (tprint>0) cout << "Pushing triangle #" << t << ", err=" << triangles[t].maxe
   << ", heap_size=" << heap_size << " now.\n"; 
 Test_Tri(elts[heap_size].tri, "after pushing message in heap::push");

 while (current>1)		// Rotate to float the largest elt.
   {
     father = current >> 1;
     if (elts[current].err() <= elts[father].err()) break;
     n_elts_pushed_onto_errs_heap_swap++;
     if (tprint>0) cout << "Heap push swapping " << current << ' ' << elts[current] << 
		       " and " << father << ' ' << elts[father] << '\n'; 
     swap(elts[current], elts[father]);
     current = father;
   }
}


// HEAP::POP

// Pop triangle with largest err, setting global vars worst_triangle and
// max_err_over_whole_TIN to the values popped off.  If the element popped is not a
// valid triangle, then repeat.  If there is no valid element, set -1 and
// 0.0 .  

// Element 0 is not used; the first heap element is 1.  This simplifies programming.
  
void Heap::pop()
{
 int c;			// Current elt's id in heap.
 int s;			// First son's id.  Second is s+1, if it exists.
 int l;			// Larger son's id.

 if (heap_size==0)
     {
	 cout << "Trying to pop an empty err heap.\n";
	 worst_triangle = -1;
	 max_err_over_whole_TIN = 0.0;
	 return ;
     }

 max_err_over_whole_TIN = elts[1].err();
 worst_triangle = elts[1].tri;
 Test_Tri(worst_triangle, "near start of heap::pop");
 n_elts_popped_from_heap++;
 errs->del(1);
 if (tprint>0) cout << "Popping triangle #" << worst_triangle << " with max_err_over_whole_TIN="
		    << max_err_over_whole_TIN << ", heap_size=" << heap_size << " now, nt=" << ntriangles << ".\n";
}


// HEAP::DEL

// Delete heap element #d, rotating if necessary.

void Heap::del(const int d)
{
    int son, current, i, l;
    
    if (d<0) 
    {  if (tprint>0) cout << "Not deleting err heap element #" << d << endl; 
       return;
    }

    i = elts[d].tri;
    swap(elts[d], elts[heap_size]);	// Move the last elt to replace the first
    heap_size--;
    n_elts_deleted_from_heap++;

    if (tprint>0) 
      cout << "Deleting heap element #" << d << ": triangle #" << i 
           << ", giving heap_size=" <<  heap_size 
           << endl;
    
    current = d;
    for (;;)		// If the first is not larger than its sons, rotate.
	{
	    son = current << 1;
	    if (son>heap_size) break;	// Don't rotate off the end of the heap.
	    if (son==heap_size)	// Is there only one son?
		{
		    if (elts[current].err() >= elts[son].err()) break;
		    l = son;
		} else {
		    if (elts[current].err() >= elts[son].err() && elts[current].err() >= elts[son+1].err()) break;
		    l = son;
		    if (elts[son+1].err() > elts[son].err()) l=son+1;
		}
	    n_swaps_while_popping_heap++;
	    if (tprint>0) cout << "Heap pop swapping " << current << ' ' << elts[current] << 
			      " and " << l << ' ' << elts[l] << '\n';
	    swap(elts[current], elts[l]);
	    current = l;
	}
}





// TEST_TRI: Test whether this is the triangle that we're watching.
// This gives us a chance to break here.

void Test_Tri(const int t, const char *msg)
{
    static int tt= -1;

    return;
    if (t == tt)
	{
	    cout << "******* I see triangle #" << t << ":" << msg << '\n';
	}
}



// TEST_PT: Test whether this is the point that we're watching.
// This gives us a chance to break here.

void Test_Pt(const int t, const char *msg)
{
    static int tt= -1;

    if (t == tt)
	{
	    cout << "******* I see point #" << t << ":" << msg << '\n';
	}
    return;
}


// WRITE HEAP_ELT

ostream& operator<<(ostream &s, Heap_elt e)
{
    return s << e.tri;
}



// LINKED_LIST CLASS ROUTINES

  void Linked_list::push(int point)		// Prefix a new point to this list
  {
    Pairindex i;

    Test_Pt(point, "Pushing onto list.");
    i = freelist.head;
    if (i==eolist) perror2("Freelist exhausted.");
    freelist.head = dotpairs[i].next;
    dotpairs[i].ptid = point;
    dotpairs[i].next = head;
    head = i;
    this->check();
  }

  int Linked_list::pop()				// Pop the first point off the list.
  {
    int i,point;
    if (head==eolist) perror2("Trying to pop a null list.");
    point = dotpairs[head].ptid;
    Test_Pt(point, "Popping off list.");
    i = dotpairs[head].next;
    dotpairs[head].ptid = -9999;
    dotpairs[head].next = freelist.head;
    freelist.head = head;
    head = i;
    freelist.check();
    return point;
  }    




// AREA of a triangle in 2D.

float Area(const Triangle t)
{
    float a0, a1, b0, b1, area;

    a0 = points[t.vert[0]].x - points[t.vert[1]].x;
    a1 = points[t.vert[0]].y - points[t.vert[1]].y;
    b0 = points[t.vert[0]].x - points[t.vert[2]].x;
    b1 = points[t.vert[0]].y - points[t.vert[2]].y;

    area = (a0*b1 - a1*b0) *0.5;

    return area;
}


// COPY_TRIANGLE Copy one triangle to another, given their ids in the
// triangles array.  Also update the heap.

void Copy_Triangle(int to, int from)
{
    if (from<0 || from>=ntriangles || to<0 || to>=ntriangles) 
	perror2("Illegal input to Copy_Triangle.");
    triangles[to] = triangles[from];
    if (triangles[to].heap_loc >= 0) 
	errs->elts[triangles[to].heap_loc].tri = to;
}

