#ifndef UNWRAP2D_H_
#define UNWRAP2D_H_
// Cloned from https://github.com/geggo/phase-unwrap
// See unwrap2D.cpp for non-commercial license terms

void unwrap2D(float* wrapped_image, float* unwrapped_image, unsigned char* input_mask,
              int image_width, int image_height,
              int wrap_around_x, int wrap_around_y);

#endif /* UNWRAP2D_H_ */
