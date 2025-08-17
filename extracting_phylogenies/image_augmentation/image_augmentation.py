import numpy as np
import matplotlib.pyplot as plt
import albumentations as A
import cv2 
import os
from io import BytesIO
from PIL import Image 


def resize_augment_wrapper(image_size):
    """
    Basic transform for training the baseline. Uniform res of <image_size>x<image_size>. Breaks aspect ratio.

    Args:
        image_size (int): height and width of resulting image in px
    """
    return A.Compose(
        [
            A.Resize(image_size,image_size),
            A.Normalize(),
        ]
    )

# Basic transform for training the baseline. Uniform res of <image_size>x<image_size>. Keeps aspect ratio!
def pad_resize_augment_wrapper(image_size, normalize, to_tensor):
    """
    Basic transform for training the baseline. Uniform res of <image_size>x<image_size> by padding instead of 
    stretching/compressing. Keeps aspect ratio of image.

    Args:
        image_size (int):  height and width of resulting image in px
        normalize (bool): 
    """
    augments= [
        # make longest side the <image_size> pixels long
        A.LongestMaxSize(max_size=image_size, interpolation=1), 
        # pad the rest for square resolution, fill with white pixels 
        A.PadIfNeeded(min_height=image_size, min_width=image_size, border_mode=0, fill=(255,255,255)), 
    ]
    if normalize:
        augments.append(A.Normalize())
    if to_tensor:
        augments.append(A.ToTensorV2())
    return A.Compose(augments)

# Resize images to uniform resolution and augment them to get more samples and make model more robust 
# towards different conditions of the image. Only used after baseline is established i.e. second training round
def dropout_augment_wrapper(image_size):
    return A.Compose(
        [   
            # make longest side the <image_size> pixels long
            A.LongestMaxSize(max_size=image_size, interpolation=1), 
            # pad the rest for square resolution
            A.PadIfNeeded(min_height=image_size, min_width=image_size, border_mode=0), 
            A.Rotate(limit=5), # too much rotation may cause taxa to be outside the visible frame
            A.OneOf([
                # holes not squre, possibly overlapping, not in grid formation
                A.CoarseDropout(num_holes_range=(1, 8), hole_height_range=(0.01, 0.1),hole_width_range=(0.01, 0.1),  p=1.0),
                # square holes in grid formation, smaller, fewer holes, non-overlapping
                A.GridDropout(ratio=0.1, p=1.0),
                # more square holes in grid format, non-overlapping
                A.GridDropout(ratio=0.2, p=1.0),
                # fewer but bigger square holes, , non-overlapping
                A.GridDropout(ratio=0.1, unit_size_range=(500, 1000),p=1.0),
            ], p=0.1), # apply one type of dropout every 10th image 
            A.Normalize(), # TODO: default values or not?
        ]
    )

def get_augmented_image(augmentation, image_path):
    """
    Augments given image and returns the bytes of the image in jpg format 
    
    Args:
        augmentation (albumentation): augmentation
        image_path (str): path to image

    Raises:
        FileNotFoundError: path not found
        IsADirectoryError: path is a dir

    Returns:
        bytes: jpg image
    """
    if not os.path.exists(image_path):
        raise FileNotFoundError(f"File not found: {image_path}")
    elif os.path.isdir(image_path):
        raise IsADirectoryError(f"Expected filepath but got directory path: {image_path}")
    # read image
    image = cv2.imread(image_path)
    # apply transform, the augmentation returns a numpy array
    augmented_nparray = augmentation(image=image)['image'] 
    # get rgb of the np array
    augmented_nparray_rgb = cv2.cvtColor(augmented_nparray, cv2.COLOR_BGR2RGB)
    # encode rgb np array to jpg
    success, jpg = cv2.imencode(".jpg", augmented_nparray_rgb)
    # return jpg encoded bytes
    if success:
        jpg_bytes = jpg.tobytes()
        return jpg_bytes
    else:
        raise ValueError("jpg encoding was not successful.")

def save_augmented_image_from_tensor(augmentation, image_path, outfile):
    if not os.path.exists(image_path):
        raise FileNotFoundError(f"File not found: {image_path}")
    elif os.path.isdir(image_path):
        raise IsADirectoryError(f"Expected filepath but got directory path: {image_path}")
    # read image
    image = cv2.imread(image_path)
    # get the augmented tensor
    augmented_tensor = augmentation(image=image)['image'] # tensor with shape (C, H, W)
    # turn tensor into numpy array
    augmented_tensor_np = np.transpose(augmented_tensor.cpu().numpy(), (1,2,0))
    # turn image into RGB format
    augmented_tensor_rgb = cv2.cvtColor(augmented_tensor_np, cv2.COLOR_BGR2RGB)
    cv2.imwrite(outfile,augmented_tensor_rgb)
    
# def main():
#     ############################## TESTING ##############################
#     # set augmentation 
#     augment_tensor = pad_resize_augment_wrapper(image_size=1024, normalize=False, to_tensor=True)
#     augment = pad_resize_augment_wrapper(image_size=1024, normalize=False, to_tensor=False)
#     # filepaths
#     test = r"C:\Users\marku\Desktop\StudiumVault\Semester6\Bachelorarbeit\Code\Extracting-Phylogenies-from-Images-using-AI\datasets\openai_finetuning_randomized\data002819460081\tree002819460081.jpg"
#     taxa_22 = r"C:\Users\marku\Desktop\StudiumVault\Semester6\Bachelorarbeit\Code\Extracting-Phylogenies-from-Images-using-AI\datasets\openai_test\data204153618932\tree204153618932.jpg"
#     taxa_22_no_labels = r"C:\Users\marku\Desktop\StudiumVault\Semester6\Bachelorarbeit\Code\Extracting-Phylogenies-from-Images-using-AI\datasets\openai_test\data204059532999\tree204059532999.jpg"
    
#     # image
#     image = get_augmented_image(augment, taxa_22_no_labels)
    
#     # view image from bytes
#     im = Image.open(BytesIO(image))

#     # Display image
#     im.show()
# if __name__ == "__main__":
#     main()