#ifndef MATHLIBRARY_H
#define MATHLIBRARY_H

#include <stdint.h>
#include <stddef.h>
#include <math.h>
#include <tgmath.h>

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

/**
 * @brief		Platform independent opaque matrix struct for holding implementation-specific handles of the matrix.
 * @attention	The matrix implementation should store the coefficients in float_t precision.
 */
typedef struct {
	void* implementation_handle;
} magneto_matrix_t;

/**
 * @brief			Platform independent function type of creating a new matrix. All values of the matrix should be 0.0f when initialized.
 * @param rows		rows of the allocated matrix.
 * @param columns	columns of the allocated matrix.
 * @retval			the created matrix.
 */
typedef magneto_matrix_t* (*new_matrix_function_ptr)(
	uint32_t rows,
	uint32_t columns
);

/**
 * @brief			Platform independent function type of deleting an existing matrix.
 * @param matrix	The matrix to be deleted.
 */
typedef void (*delete_matrix_function_ptr)(
	magneto_matrix_t* matrix
);

/**
 * @brief			Platform independent function type of getting the value of a coefficient in an existing matrix.
 * @param matrix	the matrix of the coefficient.
 * @param row		the row of the coefficient.
 * @param column	the column of the coefficient.
 * @retval			the value of the coefficient get from the matrix.
 */
typedef float_t (*get_matrix_coefficient_function_ptr)(
	magneto_matrix_t*	matrix,
	uint32_t			row,
	uint32_t			column
);

/**
 * @brief			Platform independent function type of setting the value of a coefficient in an existing matrix.
 * @param matrix	the matrix of the coefficient.
 * @param row		the row of the coefficient.
 * @param column	the column of the coefficient.
 * @param value		the value of the coefficient.
 */
typedef void (*set_matrix_coefficient_function_ptr)(
	magneto_matrix_t*	matrix,
	uint32_t			row,
	uint32_t			column,
	float_t				value
);

/**
 * @brief			Platform independent function type of adding a value to the existing value of a coefficient in an existing matrix.
 * @param matrix	the matrix of the coefficient.
 * @param row		the row of the coefficient.
 * @param column	the column of the coefficient.
 * @param value		the value to be added to the coefficient.
 */
typedef void (*add_matrix_coefficient_function_ptr)(
	magneto_matrix_t*	matrix,
	uint32_t			row,
	uint32_t			column,
	float_t				value
);

/**
 * @brief			Platform independent function type of multiplying a value to the existing value of a coefficient in an existing matrix.
 * @param matrix	the matrix of the coefficient.
 * @param row		the row of the coefficient.
 * @param column	the column of the coefficient.
 * @param value		the value to be multiplied to the coefficient.
 */
typedef void (*multiply_matrix_coefficient_function_ptr)(
	magneto_matrix_t*	matrix,
	uint32_t			row,
	uint32_t			column,
	float_t				value
);

/**
 * @brief						Platform independent function type of copying values from an existing matrix to another matrix.
 * @param source_matrix			the source matrix to be copied.
 * @param destination_matrix	destination matrix the values are copied to.
 */
typedef void (*copy_matrix_function_ptr)(
	magneto_matrix_t* source_matrix,
	magneto_matrix_t* destination_matrix
);

/**
 * @brief						Platform independent function type of copying a block of matrix from an existing matrix to a block of another existing matrix.
 * @param source_matrix			the source matrix to get the sub matrix.
 * @param from_row				the start row of the block in the source matrix.
 * @param from_column			the start column of the block in the source matrix.
 * @param to_row				the start row of the block in the destination matrix.
 * @param to_column				the start column of the block in the destination matrix.
 * @param rows					the count of rows of the block.
 * @param columns				the count of columns of the block.
 * @param destination_matrix	destination matrix the block is copied to.
 */
typedef void (*copy_matrix_block_function_ptr)(
	magneto_matrix_t*	source_matrix,
	uint32_t			from_row,
	uint32_t			from_column,
	uint32_t			to_row,
	uint32_t			to_column,
	uint32_t			rows,
	uint32_t			columns,
	magneto_matrix_t*	destination_matrix
);

/**
 * @brief						Platform independent function type of multiplying two existing matrices. (destination_matrix = left_matrix * right_matrix)
 * @param left_matrix			the left matrix of the multiplication.
 * @param right_matrix			the right matrix of the multiplication.
 * @param destination_matrix	destination matrix to hold the result matrix.
 */
typedef void (*multiply_matrix_function_ptr)(
	magneto_matrix_t* left_matrix,
	magneto_matrix_t* right_matrix,
	magneto_matrix_t* destination_matrix
);

/**
 * @brief						Platform independent function type of subtracting two existing matrices. (destination_matrix = left_matrix - right_matrix)
 * @param left_matrix			the left matrix of the subtraction.
 * @param right_matrix			the right matrix of the subtraction.
 * @param destination_matrix	destination matrix to hold the result matrix.
 */
typedef void (*subtract_matrix_function_ptr)(
	magneto_matrix_t* left_matrix,
	magneto_matrix_t* right_matrix,
	magneto_matrix_t* destination_matrix
);

/**
 * @brief				Platform independent function type of inverting an existing matrix in place.
 * @param source_matrix	the matrix to be inverted.
 */
typedef void (*invert_matrix_in_place_function_ptr)(
	magneto_matrix_t* source_matrix
);

/**
 * @brief				Platform independent function type of transposing an existing matrix in place.
 * @param source_matrix	the matrix to be transposed.
 */
typedef void (*transpose_matrix_in_place_function_ptr)(
	magneto_matrix_t* source_matrix
);

/**
 * @brief				Platform independent function type of normalizing all column vectors of an existing matrix in place.
 * @param source_matrix	the matrix to be normalized.
 */
typedef void (*normalize_matrix_in_place_function_ptr)(
	magneto_matrix_t* source_matrix
);

/**
 * @brief				Platform independent function type of setting values of all coefficients of an existing matrix to 0 in place.
 * @param source_matrix	the matrix to be set to zeros.
 */
typedef void (*set_matrix_zeros_function_ptr)(
	magneto_matrix_t* source_matrix
);

/**
 * @brief				Platform independent function type of multiplying an existing matrix with a scalar in place.
 * @param source_matrix	the matrix to be multiplied with scalar.
 * @param value			the scalar to be multiplied to the matrix.
 */
typedef void (*multiply_matrix_scalar_in_place_function_ptr)(
	magneto_matrix_t*	source_matrix,
	float_t				value
);

/**
 * @brief								Platform independent function type of solving eigenvectors and eigenvalues of an existing square matrix.
 * @attention							destination_eigenvalues_real and destination_eigenvalues_imag are column vectors (n rows 1 col matrices).
 * @param source_matrix					the source square matrix to be solved.
 * @param destination_eigenvectors_real	the real part of the eigenvectors ass column vectors packed in the matrix. The i-th column vector is the real part of the i-th eigenvector.
 * @param destination_eigenvectors_imag	the imaginary part of the eigenvectors ass column vectors packed in the matrix. The i-th column vector is the imaginary part of the i-th eigenvector.
 * @param destination_eigenvalues_real	the real part of the eigenvalues packed in the column vector. The i-th value is the real part of the i-th eigenvalue.
 * @param destination_eigenvalues_imag	the imaginary part of the eigenvalues packed in the column vector. The i-th value is the imaginary part of the i-th eigenvalue.
 * @retval								the status of the eigen solving result. ("0" = successful, other value = failed).
 */
typedef uint8_t (*solve_matrix_eigen_function_ptr)(
	magneto_matrix_t* source_matrix,
	magneto_matrix_t* destination_eigenvectors_real,
	magneto_matrix_t* destination_eigenvectors_imag,
	magneto_matrix_t* destination_eigenvalues_real,
	magneto_matrix_t* destination_eigenvalues_imag
);

/**
 * @brief Platform independent linear algebra function interface context needed by Magento.
 */
typedef struct {
	new_matrix_function_ptr							new_matrix;							/*!< Implementation function of creating a matrix. */
	delete_matrix_function_ptr						delete_matrix;						/*!< Implementation function of deleting a matrix. */
	get_matrix_coefficient_function_ptr				get_matrix_coefficient;				/*!< Implementation function of getting the value of an coefficient in a matrix. */
	set_matrix_coefficient_function_ptr				set_matrix_coefficient;				/*!< Implementation function of setting the value of an coefficient in a matrix. */
	add_matrix_coefficient_function_ptr				add_matrix_coefficient;				/*!< Implementation function of adding a value to the existing value of an coefficient in a matrix. */
	multiply_matrix_coefficient_function_ptr		multiply_matrix_coefficient;		/*!< Implementation function of multiplying a value to the existing value of an coefficient in a matrix. */
	copy_matrix_function_ptr						copy_matrix;						/*!< Implementation function of copying a matrix. */
	copy_matrix_block_function_ptr					copy_matrix_block;					/*!< Implementation function of copying a block of matrix from a matrix to block of another matrix. */
	multiply_matrix_function_ptr					multiply_matrix;					/*!< Implementation function of multiplying two matrices. */
	subtract_matrix_function_ptr					subtract_matrix;					/*!< Implementation function of subtracting two matrices. */
	invert_matrix_in_place_function_ptr				invert_matrix_in_place;				/*!< Implementation function of inverting a matrix. */
	transpose_matrix_in_place_function_ptr			transpose_matrix_in_place;			/*!< Implementation function of transposing a matrix in place. */
	normalize_matrix_in_place_function_ptr			normalize_matrix_in_place;			/*!< Implementation function of normalizing all column vectors of a matrix in place. */
	set_matrix_zeros_function_ptr					set_matrix_zeros_in_place;			/*!< Implementation function of setting all coefficients of a matrix to 0 in place. */
	multiply_matrix_scalar_in_place_function_ptr	multiply_matrix_scalar_in_place;	/*!< Implementation function of multiplying a matrix with a scalar */
	solve_matrix_eigen_function_ptr					solve_matrix_eigen;					/*!< Implementation function of solving eigenvectors and eigenvalues of a square matrix. */
} magneto_linear_algebra_context_t;

/**
 * @brief Struct for recording incoming 3-dimensional magnetometer measurement samples.
 */
typedef struct {
	magneto_matrix_t*	sample_ata_matrix;	/*!< The ATA matrix of the samples. */
	uint32_t			sample_norm_count;	/*!< The count of the samples. */
	float_t				sample_norm_sum;	/*!< The sum of length (strength) of the samples. */
} magneto_sample_container_t;

/**
 * @brief			Create a new sample container struct for recording incoming 3-dimensional magnetometer measurement samples.
 * @param context	implementation of the platform independent linear algebra function interfaces.
 * @retval			the created sample container struct.
 */
magneto_sample_container_t* magneto_new_sample_container(const magneto_linear_algebra_context_t* context);

/**
 * @brief					Reset the given sample container to clear the samples for next sample collecting.
 * @param context			implementation of the platform independent linear algebra function interfaces.
 * @param sample_container	the sample container to reset.
 */
void magneto_reset_sample_container(const magneto_linear_algebra_context_t* context, magneto_sample_container_t* sample_container);

/**
 * @brief					Delete the sample container struct.
 * @param context			implementation of the platform independent linear algebra function interfaces.
 * @param sample_container	the sample container to be deleted.
 */
void magneto_delete_sample_container(const magneto_linear_algebra_context_t* context, magneto_sample_container_t* sample_container);

/**
 * @brief					Record a sample to the struct container.
 * @param context			implementation of the platform independent linear algebra function interfaces.
 * @param sample_container	the sample container struct to be filled.
 * @param sample_x			the X-axis value of the sample.
 * @param sample_y			the Y-axis value of the sample.
 * @param sample_z			the Z-axis value of the sample.
 */
void magneto_sample(
	const	magneto_linear_algebra_context_t*	context,
			magneto_sample_container_t*			sample_container,
			float_t								sample_x,
			float_t								sample_y,
			float_t								sample_z
);

/**
 * @brief					Calculate the soft-iron calibration matrix and the hard-iron calibration vector based on the ATA matrix with samples.
 * @param context			implementation of the platform independent linear algebra function interfaces.
 * @param sample_container	the sample container struct with enough samples filled in.
 * @param soft_iron_matrix	the output matrix of the soft-iron calibration matrix, should be 3 rows and 3 cols.
 * @param hard_iron_vector	the output vector of the hard-iron calibration vector, should be 3 rows and 1 col.
 * @param reference_length	The output scalar of the reference length of the calibrated magnetometer output.
 * @retval					the status of the calculation. ("0" = no error occurred; other value = error occurs)
 */
int32_t magneto_calculate(
	const	magneto_linear_algebra_context_t*	context,
	const	magneto_sample_container_t*			sample_container,
			magneto_matrix_t*					soft_iron_matrix,
			magneto_matrix_t*					hard_iron_vector,
			float_t*							reference_length
);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // MATHLIBRARY_H
