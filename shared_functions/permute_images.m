function [permuted_images] = permute_images(images_in)

permuted_images = permute(images_in,[3,1,2]);