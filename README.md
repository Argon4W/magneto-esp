# Magneto

Magnetometer calibration using Li's ellipsoid specific fitting algorithm written in C language. This calibration is derived from the original C implementation of the [Magneto](https://sites.google.com/view/sailboatinstruments1/g-c-language-implementation). The original blog linking to the C implementation can be found [here](https://sailboatinstruments.blogspot.com/2011/09/improved-magnetometer-calibration-part.html).

## Advantage

This library separates the implementation details of linear algebra from the calibration algorithm, making it possible to attach any linear algebra library (e.g. [Eigen](https://libeigen.gitlab.io/)) to this algorithm to integrate with your project.

---

## Usage

The following example shows how to use this library to calibrate raw magnetometer output.

```c++
#include <stddef.h>
#include <math.h>

#include "magneto.h"
#include "sdkconfig.h"
#include "freertos/FreeRTOS.h"
#include "freertos/task.h"

void calibrate() {
    // Setup the math context.
    magneto_matrix_math_context_t context = {
        .new_matrix_function = ...,
        .new_matrix_array_function = ...,
        .delete_matrix_function = ...,
        .delete_matrix_array_function = ...,
        .get_matrix_element_function = ...,
        .set_matrix_element_function = ...,
        .add_matrix_element_function = ...,
        .multiply_matrix_element_function = ...,
        .copy_matrix_function = ...,
        .copy_matrix_block_function = ...,
        .copy_matrix_block2_function = ...,
        .multiply_matrix_function = ...,
        .multiply_matrix_scalar_function = ...,
        .subtract_matrix_function = ...,
        .invert_matrix_in_place_function = ...,
        .transpose_in_place_function = ...,
        .normalize_cols_in_place_function = ...,
        .negate_cols_in_place_function = ...,
        .solve_matrix_eigen_function = ...
    };
    
    // Create a sample container.
    magneto_sample_container_t* sample_container = magneto_new_sample_container(&context);
    
    for (int32_t i = 0; i < ...; i ++) {
        // Get samples from the magnetometer.
        float_t sample_x = ...;
        float_t sample_y = ...;
        float_t sample_z = ...;
        
        // Add the sample to the sample container.
        magneto_add_sample(
            &context,
            sample_container,
            sample_x,
            sample_y,
            sample_z
        );
    
        vTaskDelay(pdMS_TO_TICKS(5));
    }
    
    // Reerve space for the output.
    magneto_matrix_t* soft_iron_matrix = context.new_matrix_function(3, 3);
    magneto_matrix_t* hard_iron_vector = context.new_matrix_function(3, 1);
    
    // Calibrate.
    magneto_calculate(
        &context,
        sample_container,
        soft_iron_matrix,
        hard_iron_vector
    );
    
    // Delete the sample container.
    magneto_delete_sample_container(&context, sample_container);
    
    // Handle the calibration output like:
    float_t mag_x = ...;
    float_t mag_y = ...;
    float_t mag_z = ...;
    
    // Depends on the implementation.
    float** soft_iron_handle = soft_iron_matrix->implementation_handle;
    float** hard_iron_handle = hard_iron_vector->implementation_handle;
    
    // Apply the hard iron offset first.
    mag_x -= hard_iron_handle[0][0];
    mag_y -= hard_iron_handle[1][0];
    mag_z -= hard_iron_handle[2][0];
    
    // Then apply the soft iron matrix.
    float_t calibrated_mag_x = soft_iron_handle[0][0] * mag_x + soft_iron_handle[0][1] * mag_y + soft_iron_handle[0][2] * mag_z;
    float_t calibrated_mag_y = soft_iron_handle[1][0] * mag_x + soft_iron_handle[1][1] * mag_y + soft_iron_handle[1][2] * mag_z;
    float_t calibrated_mag_z = soft_iron_handle[2][0] * mag_x + soft_iron_handle[2][1] * mag_y + soft_iron_handle[2][2] * mag_z;
}
```