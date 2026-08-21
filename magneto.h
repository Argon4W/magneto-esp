#ifndef MATHLIBRARY_H
#define MATHLIBRARY_H

#include <string.h>
#include <stdint.h>
#include <stddef.h>
#include <math.h>
#include <tgmath.h>

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

/**
 * @brief		Platform independent opaque matrix struct for holding implementation-specific handles of the matrix.
 * @attention	The matrix implementation should store the elements in float_t precision.
 */
typedef struct {
	void* implementation_handle; /*!< Opaque handle of the matrix used by the implementation. */
} magneto_matrix_t;

/**
 * @brief	Platform independent function type of creating a new matrix. All values of the matrix should be 0.0f when initialized.
 * @retval	The created matrix.
 */
typedef magneto_matrix_t* (*new_matrix_function_ptr)(
	int32_t rows	/*!< rows of the allocated matrix. */,
	int32_t	columns	/*!< columns of the allocated matrix. */
);

/**
 * @brief	Platform independent function type of creating a new matrix array.
 * @retval	The created matrix array.
 */
typedef magneto_matrix_t* (*new_matrix_array_function_ptr)(
	int32_t count,	/*!< Count of the matrices in the array. */
	int32_t rows,	/*!< Rows of the allocated matrix in the array. */
	int32_t	columns	/*!< Columns of the allocated matrix in the array. */
);

/**
 * @brief Platform independent function type of deleting an existing matrix.
 */
typedef void (*delete_matrix_function_ptr)(
	magneto_matrix_t* matrix /*!< The matrix to be deleted. */
);

/**
 * @brief Platform independent function type of deleting an existing matrix array.
 */
typedef void (*delete_matrix_array_function_ptr)(
	magneto_matrix_t*	matrix_array,		/*!< The matrix array to be deleted. */
	int32_t				matrix_array_length	/*!< The length of the matrix array. */
);

/**
 * @brief	Platform independent function type of getting the value of an element in an existing matrix.
 * @retval	The value of the element get from the matrix.
 */
typedef float_t (*get_matrix_element_function_ptr)(
	magneto_matrix_t*	matrix,	/*!< The matrix of the element. */
	int32_t				row,	/*!< The row of the element. */
	int32_t				column	/*!< The column of the element. */
);

/**
 * @brief Platform independent function type of setting the value of an element in an existing matrix.
 */
typedef void (*set_matrix_element_function_ptr)(
	magneto_matrix_t*	matrix,	/*!< The matrix of the element. */
	int32_t				row,	/*!< The row of the element. */
	int32_t				column	/*!< The column of the element. */,
	float_t				value	/*!< The value of the element. */
);

/**
 * @brief Platform independent function type of adding a value to the existing value of an element in an existing matrix.
 */
typedef void (*add_matrix_element_function_ptr)(
	magneto_matrix_t*	matrix,	/*!< The matrix of the element. */
	int32_t				row,	/*!< The row of the element. */
	int32_t				column	/*!< The column of the element. */,
	float_t				value	/*!< The value to be added to the element. */
);

/**
 * @brief Platform independent function type of multiplying a value to the existing value of an element in an existing matrix.
 */
typedef void (*multiply_matrix_element_function_ptr)(
	magneto_matrix_t*	matrix,	/*!< The matrix of the element. */
	int32_t				row,	/*!< The row of the element. */
	int32_t				column	/*!< The column of the element. */,
	float_t				value	/*!< The value to be multiplied to the element. */
);

/**
 * @brief Platform independent function type of copying values from an existing matrix to another matrix.
 */
typedef void (*copy_matrix_function_ptr)(
	magneto_matrix_t*	source_matrix,		/*!< The source matrix to be copied. */
	magneto_matrix_t*	destination_matrix	/*!< Destination matrix to store the values, or NULL to allocate a new matrix. */
);

/**
 * @brief	Platform independent function type of copying a sub matrix from an existing matrix.
 * @retval	The result sub matrix from the source matrix.
 */
typedef magneto_matrix_t* (*copy_matrix_block_function_ptr)(
	magneto_matrix_t*	source_matrix,		/*!< The source matrix to get the sub matrix. */
	int32_t				from_row,			/*!< The start row of the sub matrix. */
	int32_t				from_column,		/*!< The start column of the sub matrix. */
	int32_t				rows,				/*!< Rows of the sub matrix. */
	int32_t				columns,			/*!< Columns of the sub matrix. */
	magneto_matrix_t*	destination_matrix	/*!< Destination matrix to store the sub matrix, or NULL to allocate a new matrix. */
);

/**
 * @brief	Platform independent function type of copying a sub matrix from an existing matrix (extended).
 * @retval	The result sub matrix from the source matrix.
 */
typedef magneto_matrix_t* (*copy_matrix_block2_function_ptr)(
	magneto_matrix_t*	source_matrix,		/*!< The source matrix to get the sub matrix. */
	int32_t				from_row,			/*!< The start row of the sub matrix. */
	int32_t				from_column,		/*!< The start column of the sub matrix. */
	int32_t				rows,				/*!< Rows of the sub matrix. */
	int32_t				columns,			/*!< Columns of the sub matrix. */
	magneto_matrix_t*	destination_matrix,	/*!< Destination matrix to store the sub matrix, or NULL to allocate a new matrix. */
	int32_t				to_row,				/*!< The first row in the destination matrix where the sub matrix stored. */
	int32_t				to_column			/*!< The first column in the destination matrix where the sub matrix stored. */
);

/**
 * @brief	Platform independent function type of multiplying two existing matrices. (destination_matrix = left_matrix * right_matrix)
 * @retval	The result multiplied matrix.
 */
typedef magneto_matrix_t* (*multiply_matrix_function_ptr)(
	magneto_matrix_t* left_matrix,			/*!< The left matrix of the multiplication. */
	magneto_matrix_t* right_matrix,			/*!< The right matrix of the multiplication. */
	magneto_matrix_t* destination_matrix	/*!< Destination matrix to hold the result matrix, or NULL to allocate a new matrix. */
);

/**
 * @brief	Platform independent function type of multiplying an existing matrix with a scalar.
 * @retval	The result multiplied matrix.
 */
typedef magneto_matrix_t* (*multiply_matrix_scalar_function_ptr)(
	float_t				value,				/*!< The scalar to be multiplied to the matrix. */
	magneto_matrix_t*	source_matrix,		/*!< The source_matrix of the multiplication. */
	magneto_matrix_t*	destination_matrix	/*!< Destination matrix to hold the result matrix, or NULL to allocate a new matrix. */
);

/**
 * @brief	Platform independent function type of subtracting two existing matrices. (destination_matrix = left_matrix - right_matrix)
 * @retval	The result multiplied matrix.
 */
typedef magneto_matrix_t* (*subtract_matrix_function_ptr)(
	magneto_matrix_t* left_matrix,			/*!< The left matrix of the subtraction. */
	magneto_matrix_t* right_matrix,			/*!< The right matrix of the subtraction. */
	magneto_matrix_t* destination_matrix	/*!< Destination matrix to hold the result matrix, or NULL to allocate a new matrix. */
);

/**
 * @brief	Platform independent function type of inverting an existing matrix in place.
 * @retval	The result inverted matrix.
 */
typedef magneto_matrix_t* (*invert_matrix_in_place_function_ptr)(
	magneto_matrix_t* source_matrix /*!< The matrix to be inverted. */
);

/**
 * @brief	Platform independent function type of transposing an existing matrix in place.
 * @retval	The result transposed matrix (identical to the input).
 */
typedef magneto_matrix_t* (*transpose_in_place_function_ptr)(
	magneto_matrix_t* source_matrix /*!< The matrix to be transposed. */
);

/**
 * @brief	Platform independent function type of normalizing all column vectors of an existing matrix in place.
 * @retval	The result normalized matrix (identical to the input).
 */
typedef magneto_matrix_t* (*normalize_cols_in_place_function_ptr)(
	magneto_matrix_t* source_matrix /*!< The matrix to be normalized. */
);

/**
 * @brief	Platform independent function type of negating all column vectors of an existing matrix in place.
 * @retval	The result negated matrix (identical to the input).
 */
typedef magneto_matrix_t* (*negate_cols_in_place_function_ptr)(
	magneto_matrix_t* source_matrix /*!< The matrix to be negated. */
);

/**
 *	@brief Result struct holding solved eigenvectors and eigenvalues of an existing matrix.
 *	@attention
 *
 *	Every element matrix in the eigenvectors array should be:
 *	[
 *		Xr, Xi,
 *		Yr, Yi,
 *		Zr, Zi
 *	];
 *
 *	The first column vector [Xr, Yr, Zr] of the matrix is the real part of the eigenvector.
 *	The second column vector [Xi, Yi, Zi] of the matrix is the imagine part of the eigenvector.
 *
 *	Every element matrix in the eigenvalues array should be:
 *	[
 *		Vr, Vi
 *	];
 *
 *	The first value (Vr) of the matrix is the real part of the eigenvalue.
 *	The second value (Vi) of the matrix is the imagine part of the eigenvalue.
 *
 *	The i-th element matrix of the eigenvalues is the eigenvalue corresponding to the i-th element matrix of the eigenvectors.
 */
typedef struct {
	uint8_t				status;			/*!< The status of the result. ("0" = successful, other value = failed). */
	magneto_matrix_t*	eigenvectors;	/*!< The solved eigenvectors array. The array size should be the same as the rows or cols of the source matrix. Every matrix in the array should be 2 cols and the same rows as the source matrix. The first column vector of the matrix is the real part of the eigen vector. The second column vector of the matrix is the imagine part of the eigen vector. */
	magneto_matrix_t*	eigenvalues;	/*!< The solved eigenvalues array, The array size should be the same as the rows or cols of the source matrix. Every matrix in the array should be 1 row and 2 cols. The first value of the matrix is the real part of the eigen value. The second value of the matrix is the imaging part of the eigen value. */
} magneto_matrix_eigen_solve_result_t;

/**
 * @brief	Platform independent function type of solving eigenvectors and eigenvalues of an existing square matrix.
 * @retval	The result struct holding the solved eigenvectors and eigenvalues.
 */
typedef magneto_matrix_eigen_solve_result_t* (*solve_matrix_eigen_function_ptr)(
	magneto_matrix_t*						source_matrix,		/*!< The square matrix to be solved. */
	magneto_matrix_eigen_solve_result_t*	destination_result	/*!< Destination result struct to hold the solved eigenvectors and eigenvalues, or NULL to allocate a new result struct and its matrices. */
);

/**
 * @brief Platform independent math function interface context needed by Magento.
 */
typedef struct {
	new_matrix_function_ptr						new_matrix_function;				/*!< Implementation function of creating a matrix. */
	new_matrix_array_function_ptr				new_matrix_array_function;			/*!< Implementation function of creating a matrix array. */
	delete_matrix_function_ptr					delete_matrix_function;				/*!< Implementation function of deleting a matrix. */
	delete_matrix_array_function_ptr			delete_matrix_array_function;		/*!< Implementation function of deleting a matrix array. */
	get_matrix_element_function_ptr				get_matrix_element_function;		/*!< Implementation function of getting the value of an element in a matrix. */
	set_matrix_element_function_ptr				set_matrix_element_function;		/*!< Implementation function of setting the value of an element in a matrix. */
	add_matrix_element_function_ptr				add_matrix_element_function;		/*!< Implementation function of adding a value to the existing value of an element in a matrix. */
	multiply_matrix_element_function_ptr		multiply_matrix_element_function;	/*!< Implementation function of multiplying a value to the existing value of an element in a matrix. */
	copy_matrix_function_ptr					copy_matrix_function;				/*!< Implementation function of copying a matrix. */
	copy_matrix_block_function_ptr				copy_matrix_block_function;			/*!< Implementation function of copying a sub matrix from a matrix. */
	copy_matrix_block2_function_ptr				copy_matrix_block2_function;		/*!< Implementation function of copying a sub matrix from a matrix. (extended) */
	multiply_matrix_function_ptr				multiply_matrix_function;			/*!< Implementation function of multiplying two matrices. */
	multiply_matrix_scalar_function_ptr			multiply_matrix_scalar_function;	/*!< Implementation function of multiplying a matrix with a scalar */
	subtract_matrix_function_ptr				subtract_matrix_function;			/*!< Implementation function of subtracting two matrices. */
	invert_matrix_in_place_function_ptr			invert_matrix_in_place_function;	/*!< Implementation function of inverting a matrix. */
	transpose_in_place_function_ptr				transpose_in_place_function;		/*!< Implementation function of transposing a matrix in place. */
	normalize_cols_in_place_function_ptr		normalize_cols_in_place_function;	/*!< Implementation function of normalizing all column vectors of a matrix in place. */
	negate_cols_in_place_function_ptr			negate_cols_in_place_function;		/*!< Implementation function of negating all column vectors of a matrix in place. */
	solve_matrix_eigen_function_ptr				solve_matrix_eigen_function;		/*!< Implementation function of solving eigenvectors and eigenvalues for a square matrix. */
} magneto_matrix_math_context_t;

/**
 * @brief Struct for containing incoming 3-dimensional magnetometer measurement samples.
 */
typedef struct {
	magneto_matrix_t*	sample_ata_matrix;	/*!< The ATA matrix of the samples. */
	float_t				sample_norm_sum;	/*!< The sum of length (strength) of the samples. */
	int32_t				sample_norm_count;	/*!< The count of the samples. */
} magneto_sample_container_t;

/**
 * @brief			Create a new sample container struct for containing incoming 3-dimensional magnetometer measurement samples.
 * @param context	implementation of the platform independent function interfaces, as pointer.
 * @retval			the created sample container struct.
 */
magneto_sample_container_t* magneto_new_sample_container(const magneto_matrix_math_context_t* context);

/**
 * @brief					Delete the sample container struct.
 * @param context			implementation of the platform independent function interfaces, as pointer.
 * @param sample_container	The sample container to be deleted.
 */
void magneto_delete_sample_container(const magneto_matrix_math_context_t* context, magneto_sample_container_t* sample_container);

/**
 * @brief					Add a sample to the struct container.
 * @param context			implementation of the platform independent function interfaces, as pointer.
 * @param sample_container	The sample container struct to be filled.
 * @param sample_x			The X-axis value of the sample.
 * @param sample_y			The Y-axis value of the sample.
 * @param sample_z			The Z-axis value of the sample.
 */
void magneto_add_sample(
	const magneto_matrix_math_context_t*	context,
	magneto_sample_container_t*				sample_container,
	const float_t							sample_x,
	const float_t							sample_y,
	const float_t							sample_z
);

/**
 * @brief					Calculate the soft-iron calibration matrix and the hard-iron calibration vector based on the ATA matrix with samples.
 * @param context			implementation of the platform independent function interfaces, as pointer.
 * @param sample_container	The sample container struct with enough samples filled in.
 * @param soft_iron_matrix	The output matrix to hold the soft-iron calibration matrix, should be 3 rows and 3 cols.
 * @param hard_iron_vector	The output vector to hold the hard-iron calibration vector, should be 3 rows and 1 cols.
 * @retval					The status of the calculation. ("0" = no error occurred; other value = error occurs)
 */
int32_t magneto_calculate(
	const magneto_matrix_math_context_t*	context,
	const magneto_sample_container_t*		sample_container,
	magneto_matrix_t*						soft_iron_matrix,
	magneto_matrix_t*						hard_iron_vector
);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // MATHLIBRARY_H
