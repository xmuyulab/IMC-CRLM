import argparse
import numpy as np
from PIL import Image
import os
from multiprocessing import Pool, cpu_count

def white_pixel_percentage(img, threshold=240):
    """Compute the percentage of white pixels in the image."""
    white_pixels = np.sum(img > threshold)
    total_pixels = img.size
    return (white_pixels / total_pixels) * 100

def normalizeStaining(img, saveFile=None, Io=240, alpha=1, beta=0.15, img_format="png"):
    ''' Normalize staining appearence of H&E stained images
    
    Example use:
        see test.py
        
    Input:
        I: RGB input image
        Io: (optional) transmitted light intensity
        
    Output:
        Inorm: normalized image
        H: hematoxylin image
        E: eosin image
    
    Reference: 
        A method for normalizing histology slides for quantitative analysis. M.
        Macenko et al., ISBI 2009
    '''
             
    HERef = np.array([[0.5626, 0.2159],
                      [0.7201, 0.8012],
                      [0.4062, 0.5581]])
        
    maxCRef = np.array([1.9705, 1.0308])
    
    # define height and width of image
    h, w, c = img.shape
    
    # reshape image
    img = img.reshape((-1,3))

    # calculate optical density
    OD = -np.log((img.astype(np.float64)+1)/Io)
    
    # remove transparent pixels
    ODhat = OD[~np.any(OD<beta, axis=1)]
        
    # compute eigenvectors
    eigvals, eigvecs = np.linalg.eigh(np.cov(ODhat.T))
    
    #eigvecs *= -1
    
    #project on the plane spanned by the eigenvectors corresponding to the two 
    # largest eigenvalues    
    That = ODhat.dot(eigvecs[:,1:3])
    
    phi = np.arctan2(That[:,1],That[:,0])
    
    minPhi = np.percentile(phi, alpha)
    maxPhi = np.percentile(phi, 100-alpha)
    
    vMin = eigvecs[:,1:3].dot(np.array([(np.cos(minPhi), np.sin(minPhi))]).T)
    vMax = eigvecs[:,1:3].dot(np.array([(np.cos(maxPhi), np.sin(maxPhi))]).T)
    
    # a heuristic to make the vector corresponding to hematoxylin first and the 
    # one corresponding to eosin second
    if vMin[0] > vMax[0]:
        HE = np.array((vMin[:,0], vMax[:,0])).T
    else:
        HE = np.array((vMax[:,0], vMin[:,0])).T
    
    # rows correspond to channels (RGB), columns to OD values
    Y = np.reshape(OD, (-1, 3)).T
    
    # determine concentrations of the individual stains
    C = np.linalg.lstsq(HE,Y, rcond=None)[0]
    
    # normalize stain concentrations
    maxC = np.array([np.percentile(C[0,:], 99), np.percentile(C[1,:],99)])
    tmp = np.divide(maxC,maxCRef)
    C2 = np.divide(C,tmp[:, np.newaxis])
    
    # recreate the image using reference mixing matrix
    Inorm = np.multiply(Io, np.exp(-HERef.dot(C2)))
    Inorm[Inorm>255] = 254
    Inorm = np.reshape(Inorm.T, (h, w, 3)).astype(np.uint8)  
    
    # unmix hematoxylin and eosin
    # H = np.multiply(Io, np.exp(np.expand_dims(-HERef[:,0], axis=1).dot(np.expand_dims(C2[0,:], axis=0))))
    # H[H>255] = 254
    # H = np.reshape(H.T, (h, w, 3)).astype(np.uint8)
    
    # E = np.multiply(Io, np.exp(np.expand_dims(-HERef[:,1], axis=1).dot(np.expand_dims(C2[1,:], axis=0))))
    # E[E>255] = 254
    # E = np.reshape(E.T, (h, w, 3)).astype(np.uint8)
    
    if saveFile is not None:
        # Check the percentage of white pixels
        if white_pixel_percentage(Inorm) <= 80:
            Image.fromarray(Inorm).save(saveFile+'.'+img_format)
        # Image.fromarray(H).save(saveFile+'_H.'+img_format)
        # Image.fromarray(E).save(saveFile+'_E.'+img_format)

    return None

def process_image(args):
    image_path, input_dir, output_dir, Io, alpha, beta, img_format = args
    # Get the relative path to maintain the directory structure
    relative_path = os.path.relpath(os.path.dirname(image_path), input_dir)
    
    # Construct a save path for the normalized image based on the original directory structure
    file = os.path.basename(image_path)
    file_without_ext = os.path.splitext(file)[0]
    save_path_base = os.path.join(output_dir, relative_path, file_without_ext)
    
    # Ensure the save directory exists
    os.makedirs(os.path.join(output_dir, relative_path), exist_ok=True)
    
    # Load and process image
    img = np.array(Image.open(image_path))
    normalizeStaining(img=img, saveFile=save_path_base, Io=Io, alpha=alpha, beta=beta, img_format=img_format)

def process_directory(input_dir, output_dir, Io, alpha, beta, img_format, nw):
    """Process all images in the directory and its subdirectories."""
    all_images = []
    for root, _, files in os.walk(input_dir):
        for file in files:
            if file.lower().endswith(('.png', '.jpg', '.jpeg', '.tiff', '.tif')):
                all_images.append(os.path.join(root, file))

    # Create argument tuples for each image
    args_list = [(img_path, input_dir, output_dir, Io, alpha, beta, img_format) for img_path in all_images]

    # Utilize all available CPU cores
    num_cores = min(cpu_count(),nw) 
    with Pool(num_cores) as pool:
        pool.map(process_image, args_list)

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputDir', type=str, help='Directory containing RGB images and their subdirectories')
    parser.add_argument('--outputDir', type=str, default='output', help='Directory to save the normalized images')
    parser.add_argument('--Io', type=int, default=240)
    parser.add_argument('--alpha', type=float, default=1)
    parser.add_argument('--beta', type=float, default=0.15)
    parser.add_argument('--format', type=str, default='png', help='Output image format (e.g. png, jpg, tiff)')
    parser.add_argument('--num_core', type=int, default=8, help='Number of core to process)')
    args = parser.parse_args()

    process_directory(args.inputDir, args.outputDir, args.Io, args.alpha, args.beta, args.format, args.num_core)

## Training sets
# nohup \
# python normalizeStaining.py --inputDir /mnt/data/lyx/IMC/WSI_HE/WSI_Data/Training_TIFF \
# --outputDir /mnt/data/lyx/IMC/WSI_HE/WSI_Data/Training_TIFF_norm \
# --Io 240 --alpha 1 --beta 0.15 --format png --num_core 16 \
# &

## Training WSI
# nohup \
# python normalizeStaining.py --inputDir /mnt/data/lyx/IMC/WSI_HE/WSI_Data/Training_TIFF \
# --outputDir /mnt/data/lyx/IMC/WSI_HE/WSI_Data/Training_TIFF_norm \
# --Io 240 --alpha 1 --beta 0.15 --format png --num_core 16 \
# &

## Test WSI
# nohup \
# python normalizeStaining.py --inputDir /mnt/data/lyx/IMC/WSI_HE/WSI_Testing/WSI_crop \
# --outputDir /mnt/data/lyx/IMC/WSI_HE/WSI_Testing/WSI_crop_norm \
# --Io 240 --alpha 1 --beta 0.15 --format png --num_core 16 \
# &