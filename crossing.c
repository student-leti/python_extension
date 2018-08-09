#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include </usr/include/python3.5/Python.h>
#include </usr/include/numpy/npy_math.h>
#include </usr/include/numpy/arrayobject.h>
#include <math.h>

static PyArrayObject *get_vector(PyArrayObject *start_point, 
								 PyArrayObject *end_point
								)
{
	PyArrayObject *vector=NULL;
	npy_intp dims[1] = {3};
	vector = (PyArrayObject *) PyArray_ZEROS(1,dims,NPY_DOUBLE,0);
	double *start_vector = (double *)PyArray_DATA(start_point);
	double *end_vector = (double *)PyArray_DATA(end_point);
	double *coord = (double *)PyArray_DATA(vector);
		
	coord[0]=end_vector[0]-start_vector[0];
	coord[1]=end_vector[1]-start_vector[1];
	coord[2]=end_vector[2]-start_vector[2];
	
	return vector;
}

static PyArrayObject *get_vector_p(PyArrayObject *vector_d1,
								   PyArrayObject *vector_d2
								  )
{
	PyArrayObject *vector_p = NULL;
	npy_intp dims[1] = {3};	   
	vector_p = (PyArrayObject *) PyArray_ZEROS(1,dims,NPY_DOUBLE,0);
	
	double *d1 = (double *)PyArray_DATA(vector_d1);
	double *d2 = (double *)PyArray_DATA(vector_d2);
	double *p = (double *)PyArray_DATA(vector_p);
	
	p[0] = (d1[1]*d2[2] - d1[2]*d2[1]);
	p[1] = -(d1[0]*d2[2] - d1[2]*d2[0]);
	p[2] = (d1[0]*d2[1] - d1[1]*d2[0]);
	
	return vector_p;
}


int check_if_parall(PyArrayObject *vector_d1,
					PyArrayObject *vector_d2
				   )
{
	double *d1 = (double *)PyArray_DATA(vector_d1);
	double *d2 = (double *)PyArray_DATA(vector_d2);
	
	if((d1[1]*d2[2] - d1[2]*d2[1])==0 &
	   (d1[0]*d2[2] - d1[2]*d2[0])==0 &
	   (d1[0]*d2[1] - d1[1]*d2[0])==0
	  )
	{
		return 1;
	}else{
		return 0;
	}
}

/*calculate the determinant of 3x3 matrix*/
double calc_determinant(PyArrayObject *row1,
						PyArrayObject *row2,
						PyArrayObject *row3
					   )
{
	double *r1 = (double *)PyArray_DATA(row1);
	double *r2 = (double *)PyArray_DATA(row2);
	double *r3 = (double *)PyArray_DATA(row3);
	
	double determinant = (r1[0]*r2[1]*r3[2] + r1[1]*r2[2]*r3[0] + r1[2]*r2[0]*r3[1] -
						  r1[2]*r2[1]*r3[0] - r1[1]*r2[0]*r3[2] - r1[0]*r2[2]*r3[1]);
	return determinant;
}

/*calculte distance between lines [a,b] and [c,d]*/
double calc_gap(PyArrayObject *vector_d1,
				PyArrayObject *vector_d2,
				PyArrayObject *vector_mm,
				PyArrayObject *vector_p
			   )
{
	double gap = 0,
		   determinant = 0,
		   scalar = 0;
		   
	double *p = (double *)PyArray_DATA(vector_p);
	
	determinant = calc_determinant(vector_d1, vector_d2, vector_mm);
	scalar = sqrtf(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)); 
	gap = fabs(determinant/scalar);
	
	return gap;
}


/*calculate distance from cross point to line segment [a,c] */
double get_distance(PyArrayObject *point_a,
					PyArrayObject *point_c,
					PyArrayObject *cross_point
				   )
{
	PyArrayObject *vector_ac = NULL,
				  *vector_ma = NULL,
				  *vector_i = NULL;
	
	double distance = 0,
		   scalar_product = 0,
		   mixed_product = 0;

	npy_intp dims[1] = {3};
				  
	vector_ac = (PyArrayObject *) PyArray_ZEROS(1,dims,NPY_DOUBLE,0);
	vector_i = (PyArrayObject *) PyArray_ZEROS(1,dims,NPY_DOUBLE,0);
	vector_ac = get_vector(point_a, point_c);
	vector_ma = get_vector(point_a, cross_point);
	
	double *ac = (double *)PyArray_DATA(vector_ac);
	double *ma = (double *)PyArray_DATA(vector_ma);
	double *i = (double *)PyArray_DATA(vector_i);
	
	i[0] = (ma[1]*ac[2] - ma[2]*ac[1]);
	i[1] = -(ma[0]*ac[2] - ma[2]*ac[0]);
	i[2] = (ma[0]*ac[1] - ma[1]*ac[0]);
		
	mixed_product = sqrtf(pow(i[0],2) + pow(i[1],2) + pow(i[2],2));
	scalar_product = sqrtf(pow(ac[0],2) + pow(ac[1],2) + pow(ac[2],2));
	distance = fabs(mixed_product/scalar_product);
	
	return distance;
}

PyArrayObject *calc_cross_point(PyArrayObject *vector_d1,
								PyArrayObject *vector_d2,
								PyArrayObject *vector_mm,
								PyArrayObject *vector_p,
								PyArrayObject *point_b,
								PyArrayObject *point_d
							   )
{
	PyArrayObject *cross_point = NULL,
				  *point_h1=NULL,
				  *point_h2=NULL;
	
	double det = 0,
		   s_det = 0,
		   t_det = 0;
		   
	npy_intp dims[1] = {3};
	
	point_h1 = (PyArrayObject *) PyArray_ZEROS(1,dims,NPY_DOUBLE,0);
	point_h2 = (PyArrayObject *) PyArray_ZEROS(1,dims,NPY_DOUBLE,0);
	cross_point = (PyArrayObject *) PyArray_ZEROS(1,dims,NPY_DOUBLE,0);
	
	double *cros_p = (double *)PyArray_DATA(cross_point);
	double *pb = (double *)PyArray_DATA(point_b);
	double *pd = (double *)PyArray_DATA(point_d);
	double *ph1 = (double *)PyArray_DATA(point_h1);
	double *ph2 = (double *)PyArray_DATA(point_h2);
	double *d1 = (double *)PyArray_DATA(vector_d1);
	double *d2 = (double *)PyArray_DATA(vector_d2);
	
	det = calc_determinant(vector_d1, vector_d2, vector_p);
	s_det = calc_determinant(vector_mm, vector_d2, vector_p);
	t_det = -calc_determinant(vector_d1, vector_mm, vector_p);
	
	ph1[0] = t_det/det*d1[0] + pb[0];
	ph1[1] = t_det/det*d1[1] + pb[1];
	ph1[2] = t_det/det*d1[2] + pb[2];
	
	ph2[0] = s_det/det*d2[0] + pd[0];
	ph2[1] = s_det/det*d2[1] + pd[1];
	ph2[2] = s_det/det*d2[2] + pd[2];
		
	cros_p[0] = (ph1[0]+ph2[0])/2;
	cros_p[1] = (ph1[1]+ph2[1])/2;
	cros_p[2] = (ph1[2]+ph2[2])/2;
	
	return cross_point;
}


static PyObject *crossing(PyObject* self, PyObject* args)
{
	PyArrayObject *point_a = NULL,
				  *point_b = NULL,
				  *point_c = NULL,
				  *point_d = NULL,
				  *vector_d1 = NULL,
				  *vector_d2 = NULL,
				  *vector_mm = NULL,
				  *vector_p = NULL,
				  *cross_point = NULL;
				  
	int flag_a = 0,
		flag_b = 0;
	
	PyArg_ParseTuple(args, "(O!O!)(O!O!)pp",
					 &PyArray_Type, &point_a,
					 &PyArray_Type, &point_b,
					 &PyArray_Type, &point_c,
					 &PyArray_Type, &point_d,
					 &flag_a,
					 &flag_b);			  
	
		
	double gap = 0,
		   determinant = 0;	

	/*calculating a direction vectors*/
	vector_d1 = get_vector(point_a, point_b);
	vector_d2 = get_vector(point_c, point_d);
	vector_mm = get_vector(point_b, point_d);
	
	/*checking if parall*/
	if(check_if_parall(vector_d1, vector_d2))
	{
		return Py_BuildValue("s","vectors are parallel");
	}
	
	/*checking if coplanar*/
	determinant = calc_determinant(vector_d1, vector_d2, vector_mm);
	if(!determinant)
	{
		return Py_BuildValue("s","vectors are not coplanar");
	}
	
	/*getin orthogonal vector*/
	vector_p = get_vector_p(vector_d1, vector_d2);
	
	/*distance between lines*/
	gap = calc_gap(vector_d1, vector_d2, vector_mm, vector_p);
	
	/*calculting cross-point*/
	cross_point = calc_cross_point(vector_d1, vector_d2, vector_mm,  vector_p, point_b, point_d);
	
	if(flag_a){}
	
	if(flag_b)
	{
		double distance=0;
		
		distance = get_distance(point_a, point_c, cross_point);
		
		return Py_BuildValue("Nfsf", cross_point, gap, flag_a ? "True" : "None", distance);
	}
	
	return Py_BuildValue("Nfss", cross_point, gap, flag_a ? "True" : "None", "None");
}

static PyMethodDef module_methods[] = {
	{"crossing", crossing,
		METH_VARARGS, "summing arrays"},
		{NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
	PyModuleDef_HEAD_INIT,
	"summing_numpy", /*name of module*/
	"summing arrays",
	-1,
	module_methods
};

PyMODINIT_FUNC PyInit_summing_numpy(void){
	PyObject *m;
	m = PyModule_Create(&moduledef);
	/*!! Import NUMPY*/
	import_array();
	
	return m;
}
