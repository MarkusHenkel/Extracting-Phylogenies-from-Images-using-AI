import albumentations as A
import cv2 
import os
from io import BytesIO
from PIL import Image 


def resize_augment_wrapper(height=512,width=512, normalize=False, to_tensor=False):
    """
    Basic transform for training the baseline. Uniform res of <image_size>x<image_size>. Breaks aspect ratio.

    Args:
        image_size (int): height and width of resulting image in px
        normalize (bool): apply default normalization of pixel values 
        to_tensor (bool): Converts images/masks to PyTorch Tensors, inheriting from BasicTransform. For images: 
        If input is in HWC format, converts to PyTorch CHW format. If input is in HW format, converts to PyTorch 1HW format 
        (adds channel dimension)
    """
    augments= [A.Resize(height,width)]
    if normalize:
        augments.append(A.Normalize())
    if to_tensor:
        augments.append(A.ToTensorV2())
    return A.Compose(augments)

# Basic transform for training the baseline. Uniform res of <image_size>x<image_size>. Keeps aspect ratio!
def pad_resize_augment_wrapper(image_size=512, normalize=False, to_tensor=False):
    """
    Basic transform for training the baseline. Uniform res of <image_size>x<image_size> by padding instead of 
    stretching/compressing. Keeps aspect ratio of image.

    Args:
        image_size (int):  height and width of resulting image in px
        normalize (bool): apply default normalization of pixel values 
        to_tensor (bool): Converts images/masks to PyTorch Tensors, inheriting from BasicTransform. For images: 
        If input is in HWC format, converts to PyTorch CHW format. If input is in HW format, converts to PyTorch 1HW format 
        (adds channel dimension)
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
def dropout_augment_wrapper(image_size=512, dropout_probability=0.1, normalize=False, to_tensor=False):
    augments = [
        # make longest side the <image_size> pixels long
        A.LongestMaxSize(max_size=image_size, interpolation=1), 
        # pad the rest for square resolution
        A.PadIfNeeded(min_height=image_size, min_width=image_size, border_mode=0), 
        A.Rotate(limit=5), # too much rotation may cause taxa to be outside the visible frame
        A.OneOf([
            # holes not square, possibly overlapping, not in grid formation
            A.CoarseDropout(num_holes_range=(1, 8), hole_height_range=(0.01, 0.1),hole_width_range=(0.01, 0.1),  p=1.0),
            # square holes in grid formation, smaller, fewer holes, non-overlapping
            A.GridDropout(ratio=0.05, p=1.0),
            # more square holes in grid format, non-overlapping
            A.GridDropout(ratio=0.1, p=1.0),
            # fewer but bigger square holes, non-overlapping
            A.GridDropout(ratio=0.05, unit_size_range=(200, 500),p=1.0),
        ], p=dropout_probability), # apply one type of dropout every 10th image 
    ]
    if normalize:
        augments.append(A.Normalize())
    if to_tensor:
        augments.append(A.ToTensorV2())
    return A.Compose(augments)

def get_augmented_image(augmentation, image_path):
    """
    Augments given image and returns a PIL object of the image
    
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
        image = Image.open(BytesIO(jpg_bytes))
        return image
    else:
        raise ValueError("jpg encoding was not successful.")
    
def show_augmented_image(augmentation, image_path):
    image = get_augmented_image(augmentation=augmentation,image_path=image_path)
    image.show()
    
def save_augmented_image(augmentation, image_path, outfile):
    image = get_augmented_image(augmentation=augmentation, image_path=image_path)
    image.save(outfile)