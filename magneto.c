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

magneto_sample_container_t* magneto_new_sample_container(const magneto_matrix_math_context_t* context) {
	magneto_sample_container_t* sample_container = (magneto_sample_container_t*) malloc(sizeof(magneto_sample_container_t));

	sample_container->sample_ata_matrix	= context->new_matrix_function(10, 10);
	sample_container->sample_norm_sum	= 0.0f;
	sample_container->sample_norm_count	= 0;

	return sample_container;
}

void magneto_delete_sample_container(const magneto_matrix_math_context_t* context, magneto_sample_container_t* sample_container) {
	context->delete_matrix_function(sample_container->sample_ata_matrix);
	free(sample_container);
}

void magneto_add_sample(
	const magneto_matrix_math_context_t*	context,
	magneto_sample_container_t*				sample_container,
	const float_t							sample_x,
	const float_t							sample_y,
	const float_t							sample_z
) {
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

	sample_container->sample_norm_count	++;
	sample_container->sample_norm_sum	+= sqrt(
		sample_x * sample_x +
		sample_y * sample_y +
		sample_z * sample_z
	);

	for		(int32_t row = 0; row < 10; row ++) {
		for	(int32_t col = 0; col < 10; col ++) {
			context->add_matrix_element_function(
				sample_container->sample_ata_matrix,
				row,
				col,
				d[row] * d[col]
			);
		}
	}
}

int32_t magneto_calculate(
	const magneto_matrix_math_context_t*	context,
	const magneto_sample_container_t*		sample_container,
	magneto_matrix_t*						soft_iron_matrix,
	magneto_matrix_t*						hard_iron_vector
) {
	if (sample_container->sample_norm_count == 0) {
		return -1;
	}

	magneto_matrix_t* S11	= context->new_matrix_function(6, 6);
	magneto_matrix_t* S12	= context->new_matrix_function(6, 4);
	magneto_matrix_t* S12t	= context->new_matrix_function(4, 6);
	magneto_matrix_t* S22	= context->new_matrix_function(4, 4);

	context->copy_matrix_block_function(sample_container->sample_ata_matrix, 0, 0, 6, 6, S11);
	context->copy_matrix_block_function(sample_container->sample_ata_matrix, 0, 6, 6, 4, S12);
	context->copy_matrix_block_function(sample_container->sample_ata_matrix, 6, 0, 4, 6, S12t);
	context->copy_matrix_block_function(sample_container->sample_ata_matrix, 6, 6, 4, 4, S22);

	context->invert_matrix_in_place_function(S22);

	magneto_matrix_t* S22a = context->new_matrix_function(4, 6);
	magneto_matrix_t* S22b = context->new_matrix_function(6, 6);

	context->multiply_matrix_function	(S22, S12t, S22a); // Calculate S22a = S22_1 * S12t; 4*6 = 4x4 * 4x6; C = AB
	context->multiply_matrix_function	(S12, S22a, S22b); // Then calculate S22b = S12 * S22a ( 6x6 = 6x4 * 4x6)
	context->delete_matrix_function		(S22);
	context->delete_matrix_function		(S12t);
	context->delete_matrix_function		(S12);

	// Calculate SS = S11 - S22b.
	magneto_matrix_t* SS	= context->new_matrix_function(6, 6);
	magneto_matrix_t* E		= context->new_matrix_function(6, 6);
	magneto_matrix_t* C		= context->new_matrix_function(6, 6);

	// Fill the matrix C.
	for		(int32_t row = 0; row < 6; row ++) {
		for	(int32_t col = 0; col < 6; col ++) {
			context->set_matrix_element_function(C, row, col, C_value[row * 6 + col]);
		}
	}

	context->subtract_matrix_function	(S11,	S22b,	SS);
	context->multiply_matrix_function	(C,		SS,		E);
	context->delete_matrix_function		(S22b);
	context->delete_matrix_function		(S11);
	context->delete_matrix_function		(SS);
	context->delete_matrix_function		(C);

	magneto_matrix_t* eigenvectors	= context->new_matrix_array_function(6, 6, 2);
	magneto_matrix_t* eigenvalues	= context->new_matrix_array_function(6, 1, 2);

	magneto_matrix_eigen_solve_result_t eigen_solve_result = {
		.status			= 0,
		.eigenvectors	= eigenvectors,
		.eigenvalues	= eigenvalues
	};

	context->solve_matrix_eigen_function(E, &eigen_solve_result);
	context->delete_matrix_function		(E);

	if (eigen_solve_result.status != 0) {
		context->delete_matrix_function			(S22a);
		context->delete_matrix_array_function	(eigenvectors,	6);
		context->delete_matrix_array_function	(eigenvalues,	6);

		return eigen_solve_result.status;
	}

	int32_t				max_index = 0;
	magneto_matrix_t*	max_value = &eigenvalues[0];

	for (int32_t i = 1; i < 6; i++) {
		magneto_matrix_t*	new_value	= &eigenvalues[i];
		const float_t		new_real	= context->get_matrix_element_function(&eigenvalues[i],	0, 0);
		const float_t		max_real	= context->get_matrix_element_function(max_value,		0, 0);

		if(new_real > max_real) {
			max_value = new_value;
			max_index = i;
		}
	}

	magneto_matrix_t* v1 = context->new_matrix_function(6, 1);
	magneto_matrix_t* v2 = context->new_matrix_function(4, 1);

	context->copy_matrix_block_function(
		/* source_matrix		= */ &eigenvectors[max_index],
		/* from_row				= */ 0,
		/* from_col				= */ 0,
		/* rows					= */ 6,
		/* cols					= */ 1,
		/* destination_matrix	= */ v1
	);

	context->delete_matrix_array_function(eigenvectors,	6);
	context->delete_matrix_array_function(eigenvalues,	6);

	context->normalize_cols_in_place_function(v1);

	if (context->get_matrix_element_function(v1, 0, 0) < 0.0f) {
		context->negate_cols_in_place_function(v1);
	}

	// Calculate v2 = S22a * v1 ( 4x1 = 4x6 * 6x1).
	context->multiply_matrix_function		(S22a, v1, v2);
	context->delete_matrix_function			(S22a);
	context->negate_cols_in_place_function	(v2);

	magneto_matrix_t* Q		= context->new_matrix_function(3, 3);
	magneto_matrix_t* Qi	= context->new_matrix_function(3, 3);
	magneto_matrix_t* U		= context->new_matrix_function(3, 1);

	context->set_matrix_element_function(Q, 0, 0, context->get_matrix_element_function(v1, 0, 0));
	context->set_matrix_element_function(Q, 0, 1, context->get_matrix_element_function(v1, 5, 0));
	context->set_matrix_element_function(Q, 0, 2, context->get_matrix_element_function(v1, 4, 0));

	context->set_matrix_element_function(Q, 1, 0, context->get_matrix_element_function(v1, 5, 0));
	context->set_matrix_element_function(Q, 1, 1, context->get_matrix_element_function(v1, 1, 0));
	context->set_matrix_element_function(Q, 1, 2, context->get_matrix_element_function(v1, 3, 0));

	context->set_matrix_element_function(Q, 2, 0, context->get_matrix_element_function(v1, 4, 0));
	context->set_matrix_element_function(Q, 2, 1, context->get_matrix_element_function(v1, 3, 0));
	context->set_matrix_element_function(Q, 2, 2, context->get_matrix_element_function(v1, 2, 0));

	context->copy_matrix_block_function(
		/* source_matrix		= */ v2,
		/* from_row				= */ 0,
		/* from_col				= */ 0,
		/* rows					= */ 3,
		/* cols					= */ 1,
		/* destination_matrix	= */ U
	);

	const float_t J = context->get_matrix_element_function(v2, 3, 0);

	context->delete_matrix_function(v1);
	context->delete_matrix_function(v2);

	context->copy_matrix_function			(Q, Qi);
	context->invert_matrix_in_place_function(Qi);

	// Calculate B = Q-1 * U ( 3x1 = 3x3 * 3x1) (B = hard_iron_vector).
	context->multiply_matrix_function		(Qi, U, hard_iron_vector);
	context->delete_matrix_function			(Qi);
	context->delete_matrix_function			(U);
	context->negate_cols_in_place_function	(hard_iron_vector);

	magneto_matrix_t* QB = context->new_matrix_function(3, 1);

	// First calculate QB = Q * B   ( 3x1 = 3x3 * 3x1).
	context->multiply_matrix_function(Q, hard_iron_vector, QB);

	// 1*1 matrix.
	magneto_matrix_t* Bt	= context->new_matrix_function(3, 1);
	magneto_matrix_t* BtQB	= context->new_matrix_function(1, 1);

	// Then calculate BtQB = BT * QB    ( 1x1 = 1x3 * 3x1).
	context->copy_matrix_function		(hard_iron_vector, Bt);
	context->transpose_in_place_function(Bt);
	context->multiply_matrix_function	(Bt, QB, BtQB);

	// Calculate hmb = sqrt(BtQb - J).
	const float_t hmb = sqrt(context->get_matrix_element_function(BtQB, 0, 0) - J);

	context->delete_matrix_function(QB);
	context->delete_matrix_function(Bt);
	context->delete_matrix_function(BtQB);

	// Calculate SQ, the square root of matrix Q.
	magneto_matrix_t* eigenvectors2	= context->new_matrix_array_function(3, 3, 2);
	magneto_matrix_t* eigenvalues2	= context->new_matrix_array_function(3, 1, 2);

	magneto_matrix_eigen_solve_result_t eigen_solve_result2 = {
		.status			= 0,
		.eigenvectors	= eigenvectors2,
		.eigenvalues	= eigenvalues2
	};

	context->solve_matrix_eigen_function(Q, &eigen_solve_result2);
	context->delete_matrix_function		(Q);

	if (eigen_solve_result2.status != 0) {
		context->delete_matrix_array_function(eigenvectors2,	3);
		context->delete_matrix_array_function(eigenvalues2,		3);

		return eigen_solve_result2.status;
	}

	// Normalize eigenvectors.
	for (int32_t i = 0; i < 3; i ++) {
		context->normalize_cols_in_place_function(&eigenvectors2[i]);
	}

	magneto_matrix_t* Dz = context->new_matrix_function(3, 3);

	context->set_matrix_element_function(Dz, 0, 0, sqrt(context->get_matrix_element_function(&eigenvalues2[0], 0, 0)));
	context->set_matrix_element_function(Dz, 1, 1, sqrt(context->get_matrix_element_function(&eigenvalues2[1], 0, 0)));
	context->set_matrix_element_function(Dz, 2, 2, sqrt(context->get_matrix_element_function(&eigenvalues2[2], 0, 0)));

	magneto_matrix_t* eigenvectors_real	= context->new_matrix_function(3, 3);
	magneto_matrix_t* vdz				= context->new_matrix_function(3, 3);
	magneto_matrix_t* SQ				= context->new_matrix_function(3, 3);

	context->copy_matrix_block2_function(&eigenvectors2[0], 0, 0, 3, 1, eigenvectors_real, 0, 0);
	context->copy_matrix_block2_function(&eigenvectors2[1], 0, 0, 3, 1, eigenvectors_real, 0, 1);
	context->copy_matrix_block2_function(&eigenvectors2[2], 0, 0, 3, 1, eigenvectors_real, 0, 2);

	context->delete_matrix_array_function(eigenvectors2,	3);
	context->delete_matrix_array_function(eigenvalues2,		3);

	const float_t hm = sample_container->sample_norm_sum / (float) sample_container->sample_norm_count;

	context->multiply_matrix_function		(eigenvectors_real, Dz, vdz);
	context->transpose_in_place_function	(eigenvectors_real);
	context->multiply_matrix_function		(vdz, eigenvectors_real, SQ);
	context->multiply_matrix_scalar_function(hm / hmb, SQ, soft_iron_matrix);

	context->delete_matrix_function(Dz);
	context->delete_matrix_function(SQ);
	context->delete_matrix_function(vdz);
	context->delete_matrix_function(eigenvectors_real);

	return 0;
}