#include <magneto.h>

// Pre-inverted constraint matrix C
const float_t C_value[6 * 6] = {
	0.0f,	0.5f,	0.5f,	0.0f,	0.0f,	0.0f,
	0.5f,	0.0f,	0.5f,	0.0f,	0.0f,	0.0f,
	0.5f,	0.5f,	0.0f,	0.0f,	0.0f,	0.0f,
	0.0f,	0.0f,	0.0f,	-0.25f,	0.0f,	0.0f,
	0.0f,	0.0f,	0.0f,	0.0f,	-0.25f,	0.0f,
	0.0f,	0.0f,	0.0f,	0.0f,	0.0f,	-0.25f
};

magneto_sample_container_t* magneto_new_sample_container(const magneto_linear_algebra_context_t* context) {
	// Allocate the sample container in heap.
	magneto_sample_container_t* sample_container = (magneto_sample_container_t*) malloc(sizeof(magneto_sample_container_t));

	// Initialize the coefficients of the sample container.
	sample_container->sample_ata_matrix	= context->new_matrix(10U, 10U);	// Create the ATA matrix (information matrix/cross-product matrx) for the fitting.
	sample_container->sample_norm_sum	= 0.0f;								// Reset the sum of magnetometer sample strengths to 0.0f to avoid undefined values.
	sample_container->sample_norm_count	= 0U;								// Reset the count of magnetometer samples to 0 to avoid undefined valueus.

	return sample_container;
}

void magneto_delete_sample_container(const magneto_linear_algebra_context_t* context, magneto_sample_container_t* sample_container) {
	// First delete the ATA matrix inside the sample container to avoid leaking.
	context->delete_matrix(sample_container->sample_ata_matrix);
	// Then we can safely delete the container itself.
	free(sample_container);
}

void magneto_sample(
	const	magneto_linear_algebra_context_t*	context,
			magneto_sample_container_t*			sample_container,
	const	float_t								sample_x,
	const	float_t								sample_y,
	const	float_t								sample_z
) {
	// Pre-calculate all coefficients of the sample.
	const float_t d[10] = {
		sample_x * sample_x,
		sample_y * sample_y,
		sample_z * sample_z,
		2.0f * sample_y * sample_z,
		2.0f * sample_x * sample_z,
		2.0f * sample_x * sample_y,
		2.0f * sample_x,
		2.0f * sample_y,
		2.0f * sample_z,
		1.0f
	};

	// Update the sample container.
	sample_container->sample_norm_count	++;			// UUpdate the sample count.
	sample_container->sample_norm_sum	+= sqrt(	// Accumulate the strengths of the samples.
		sample_x * sample_x +
		sample_y * sample_y +
		sample_z * sample_z
	);

	// Accumulate the sample to the ATA matrix of the sample container.
	for		(uint32_t row = 0; row < 10U; row ++) {
		for	(uint32_t col = 0; col < 10U; col ++) {
			// Add the coefficients to the ATA matrix.
			context->add_matrix_coefficient(
				/* matrix	= */ sample_container->sample_ata_matrix,
				/* row		= */ row,
				/* column	= */ col,
				/* value	= */ d[row] * d[col]
			);
		}
	}
}

int32_t magneto_calculate(
	const	magneto_linear_algebra_context_t*	context,
	const	magneto_sample_container_t*			sample_container,
			magneto_matrix_t*					soft_iron_matrix,
			magneto_matrix_t*					hard_iron_vector
) {
	if (sample_container->sample_norm_count == 0) {
		return -1;
	}

	magneto_matrix_t* S11	= context->new_matrix(6U, 6U);
	magneto_matrix_t* S12	= context->new_matrix(6U, 4U);
	magneto_matrix_t* S12t	= context->new_matrix(4U, 6U);
	magneto_matrix_t* S22	= context->new_matrix(4U, 4U);

	context->copy_matrix_block(sample_container->sample_ata_matrix, 0U, 0U, 0U, 0U, 6U, 6U, S11);
	context->copy_matrix_block(sample_container->sample_ata_matrix, 0U, 6U, 0U, 0U, 6U, 4U, S12);
	context->copy_matrix_block(sample_container->sample_ata_matrix, 6U, 0U, 0U, 0U, 4U, 6U, S12t);
	context->copy_matrix_block(sample_container->sample_ata_matrix, 6U, 6U, 0U, 0U, 4U, 4U, S22);

	context->invert_matrix_in_place(S22);

	magneto_matrix_t* S22a = context->new_matrix(4U, 6U);
	magneto_matrix_t* S22b = context->new_matrix(6U, 6U);

	context->multiply_matrix(S22, S12t, S22a); // Calculate S22a = S22_1 * S12t; 4*6 = 4x4 * 4x6; C = AB
	context->multiply_matrix(S12, S22a, S22b); // Then calculate S22b = S12 * S22a ( 6x6 = 6x4 * 4x6)

	context->delete_matrix(S22);
	context->delete_matrix(S12t);
	context->delete_matrix(S12);

	// Calculate SS = S11 - S22b.
	magneto_matrix_t* SS	= context->new_matrix(6U, 6U);
	magneto_matrix_t* E		= context->new_matrix(6U, 6U);
	magneto_matrix_t* C		= context->new_matrix(6U, 6U);

	// Fill the constraint matrix C.
	for		(uint32_t row = 0U; row < 6U; row ++) {
		for	(uint32_t col = 0U; col < 6U; col ++) {
			context->set_matrix_coefficient(C, row, col, C_value[row * 6U + col]);
		}
	}

	context->subtract_matrix(S11,	S22b,	SS);
	context->multiply_matrix(C,		SS,		E);

	context->delete_matrix(S22b);
	context->delete_matrix(S11);
	context->delete_matrix(SS);
	context->delete_matrix(C);

	magneto_matrix_t* eigenvectors_real	= context->new_matrix(6U, 6U);
	magneto_matrix_t* eigenvectors_imag	= context->new_matrix(6U, 6U);
	magneto_matrix_t* eigenvalues_real	= context->new_matrix(6U, 1U);
	magneto_matrix_t* eigenvalues_imag	= context->new_matrix(6U, 1U);

	const uint8_t result = context->solve_matrix_eigen(
		/* source_matrix					= */ E,
		/* destination_eigenvectors_real	= */ eigenvectors_real,
		/* destination_eigenvectors_imag	= */ eigenvectors_imag,
		/* destination_eigenvalues_real		= */ eigenvalues_real,
		/* destination_eigenvalues_imag		= */ eigenvalues_imag
	);

	context->delete_matrix(E);
	context->delete_matrix(eigenvectors_imag);
	context->delete_matrix(eigenvalues_imag);

	if (result != 0U) {
		context->delete_matrix(S22a);
		context->delete_matrix(eigenvectors_real);
		context->delete_matrix(eigenvalues_real);

		return result;
	}

	uint32_t	max_index = 0;
	float_t		max_value = context->get_matrix_coefficient(eigenvalues_real, 0U, 0U);

	for (uint32_t i = 1U; i < 6U; i++) {
		const float_t new_value = context->get_matrix_coefficient(eigenvalues_real, i, 0U);

		if (new_value > max_value) {
			max_value = new_value;
			max_index = i;
		}
	}

	magneto_matrix_t* v1 = context->new_matrix(6U, 1U);
	magneto_matrix_t* v2 = context->new_matrix(4U, 1U);

	context->copy_matrix_block(
		/* source_matrix		= */ eigenvectors_real,
		/* from_row				= */ 0U,
		/* from_col				= */ max_index,
		/* to_row				= */ 0U,
		/* to_col				= */ 0U,
		/* rows					= */ 6U,
		/* cols					= */ 1U,
		/* destination_matrix	= */ v1
	);

	context->delete_matrix(eigenvectors_real);
	context->delete_matrix(eigenvalues_real);

	context->normalize_matrix_in_place(v1);

	if (context->get_matrix_coefficient(v1, 0U, 0U) < 0.0f) {
		context->multiply_matrix_scalar_in_place(v1, -1);
	}

	// Calculate v2 = S22a * v1 ( 4x1 = 4x6 * 6x1).
	context->multiply_matrix	(S22a, v1, v2);
	context->delete_matrix		(S22a);

	context->multiply_matrix_scalar_in_place(v2, -1.0f);

	magneto_matrix_t* Q		= context->new_matrix(3U, 3U);
	magneto_matrix_t* Qi	= context->new_matrix(3U, 3U);
	magneto_matrix_t* U		= context->new_matrix(3U, 1U);

	context->set_matrix_coefficient(Q, 0U, 0U, context->get_matrix_coefficient(v1, 0U, 0U));
	context->set_matrix_coefficient(Q, 0U, 1U, context->get_matrix_coefficient(v1, 5U, 0U));
	context->set_matrix_coefficient(Q, 0U, 2U, context->get_matrix_coefficient(v1, 4U, 0U));

	context->set_matrix_coefficient(Q, 1U, 0U, context->get_matrix_coefficient(v1, 5U, 0U));
	context->set_matrix_coefficient(Q, 1U, 1U, context->get_matrix_coefficient(v1, 1U, 0U));
	context->set_matrix_coefficient(Q, 1U, 2U, context->get_matrix_coefficient(v1, 3U, 0U));

	context->set_matrix_coefficient(Q, 2U, 0U, context->get_matrix_coefficient(v1, 4U, 0U));
	context->set_matrix_coefficient(Q, 2U, 1U, context->get_matrix_coefficient(v1, 3U, 0U));
	context->set_matrix_coefficient(Q, 2U, 2U, context->get_matrix_coefficient(v1, 2U, 0U));

	context->copy_matrix_block(
		/* source_matrix		= */ v2,
		/* from_row				= */ 0U,
		/* from_col				= */ 0U,
		/* to_row				= */ 0U,
		/* to_col				= */ 0U,
		/* rows					= */ 3U,
		/* cols					= */ 1U,
		/* destination_matrix	= */ U
	);

	const float_t J = context->get_matrix_coefficient(v2, 3U, 0U);

	context->delete_matrix(v1);
	context->delete_matrix(v2);

	context->copy_matrix			(Q, Qi);
	context->invert_matrix_in_place	(Qi);

	// Calculate B = Q-1 * U ( 3x1 = 3x3 * 3x1) (B = hard_iron_vector).
	context->multiply_matrix					(Qi, U, hard_iron_vector);
	context->multiply_matrix_scalar_in_place	(hard_iron_vector, -1.0f);

	context->delete_matrix(Qi);
	context->delete_matrix(U);

	magneto_matrix_t* QB = context->new_matrix(3U, 1U);

	// First calculate QB = Q * B   ( 3x1 = 3x3 * 3x1).
	context->multiply_matrix(Q, hard_iron_vector, QB);

	// 1*1 matrix.
	magneto_matrix_t* Bt	= context->new_matrix(3U, 1U);
	magneto_matrix_t* BtQB	= context->new_matrix(1U, 1U);

	// Then calculate BtQB = BT * QB    ( 1x1 = 1x3 * 3x1).
	context->copy_matrix				(hard_iron_vector, Bt);
	context->transpose_matrix_in_place	(Bt);
	context->multiply_matrix			(Bt, QB, BtQB);

	// Calculate hmb = sqrt(BtQb - J).
	const float_t hmb = sqrt(context->get_matrix_coefficient(BtQB, 0U, 0U) - J);

	context->delete_matrix(QB);
	context->delete_matrix(Bt);
	context->delete_matrix(BtQB);

	// Calculate SQ, the square root of matrix Q.
	magneto_matrix_t* eigenvectors_real2	= context->new_matrix(3U, 3U);
	magneto_matrix_t* eigenvectors_imag2	= context->new_matrix(3U, 3U);
	magneto_matrix_t* eigenvalues_real2		= context->new_matrix(3U, 1U);
	magneto_matrix_t* eigenvalues_imag2		= context->new_matrix(3U, 1U);

	const uint8_t result2 = context->solve_matrix_eigen(
		/* source_matrix					= */ Q,
		/* destination_eigenvectors_real	= */ eigenvectors_real2,
		/* destination_eigenvectors_imag	= */ eigenvectors_imag2,
		/* destination_eigenvalues_real		= */ eigenvalues_real2,
		/* destination_eigenvalues_imag		= */ eigenvalues_imag2
	);

	context->delete_matrix(Q);
	context->delete_matrix(eigenvectors_imag2);
	context->delete_matrix(eigenvalues_imag2);

	if (result2 != 0U) {
		context->delete_matrix(eigenvectors_real2);
		context->delete_matrix(eigenvalues_real2);

		return result2;
	}

	// Normalize all column eigenvectors in the matrix.
	context->normalize_matrix_in_place(eigenvectors_real2);

	magneto_matrix_t* Dz = context->new_matrix(3U, 3U);

	context->set_matrix_coefficient(Dz, 0U, 0U, sqrt(context->get_matrix_coefficient(eigenvalues_real2, 0U, 0U)));
	context->set_matrix_coefficient(Dz, 1U, 1U, sqrt(context->get_matrix_coefficient(eigenvalues_real2, 1U, 0U)));
	context->set_matrix_coefficient(Dz, 2U, 2U, sqrt(context->get_matrix_coefficient(eigenvalues_real2, 2U, 0U)));

	magneto_matrix_t* vdz	= context->new_matrix(3U, 3U);
	magneto_matrix_t* SQ	= context->new_matrix(3U, 3U);

	const float_t hm = sample_container->sample_norm_sum / (float) sample_container->sample_norm_count;

	context->multiply_matrix			(eigenvectors_real2, Dz, vdz);
	context->transpose_matrix_in_place	(eigenvectors_real2);
	context->multiply_matrix			(vdz, eigenvectors_real2, SQ);

	context->multiply_matrix_scalar_in_place(SQ, hm / hmb);
	context->copy_matrix					(SQ, soft_iron_matrix);

	context->delete_matrix(Dz);
	context->delete_matrix(SQ);
	context->delete_matrix(vdz);
	context->delete_matrix(eigenvectors_real2);
	context->delete_matrix(eigenvalues_real2);

	return 0;
}