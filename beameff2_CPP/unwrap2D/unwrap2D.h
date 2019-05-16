#ifndef UNWRAP2D_H_
#define UNWRAP2D_H_
// cloned from https://github.com/scikit-image/scikit-image/blob/master/skimage/restoration/unwrap_2d_ljmu.c
// having license (Modified BSD) from https://github.com/scikit-image/scikit-image

void unwrap2D(double *wrapped_image, double *UnwrappedImage,
              unsigned char *input_mask, int image_width, int image_height,
              int wrap_around_x, int wrap_around_y,
              char use_seed = 0, unsigned int seed = 0);

#endif /* UNWRAP2D_H_ */
