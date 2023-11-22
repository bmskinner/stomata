import os, json, cv2, numpy as np
import matplotlib
import matplotlib.pyplot as plt

import torch
from torch.utils.data import Dataset, DataLoader

import torchvision
from torchvision.models.detection.rpn import AnchorGenerator
from torchvision.transforms import functional as F

import albumentations as A # Library for augmentations

# https://github.com/pytorch/vision/tree/main/references/detection
import transforms, utils, engine, train
from utils import collate_fn
from engine import train_one_epoch, evaluate

matplotlib.use('TkAgg') # ensure plot window can spawn

MODEL_NAME = 'keypointsrcnn_weights_3.pth'
data_path = './stomatal_images_trainable'
KEYPOINTS_FOLDER_TRAIN = os.path.join(data_path, 'train')
KEYPOINTS_FOLDER_TEST = os.path.join(data_path, 'test')

# Define the function for augmentations for training. 
# This function will apply different transforms to the images before each training iteration. 
#Among such transforms there could be a random change of brightness and contrast or an 
# image rotation by 90 degrees random number of times.

#Thus, we essentially “create new images”, which are in some ways differ from 
# original ones, but still perfectly suitable for training our model.

def train_transform():
    return A.Compose([
        A.Sequential([
            A.RandomRotate90(p=1), # Random rotation of an image by 90 degrees zero or more times
            A.RandomBrightnessContrast(brightness_limit=0.1, contrast_limit=0.1, 
            brightness_by_max=True, always_apply=False, p=1), # Random change of brightness & contrast
        ], p=1)
    ],
    keypoint_params=A.KeypointParams(format='xy'), # More about keypoint formats used in albumentations library read at https://albumentations.ai/docs/getting_started/keypoints_augmentation/
    bbox_params=A.BboxParams(format='pascal_voc', label_fields=['bboxes_labels']) # Bboxes should have labels, read more at https://albumentations.ai/docs/getting_started/bounding_boxes_augmentation/
    )


# The dataset should inherit from the standard torch.utils.data.Dataset class, 
# and __getitem__ should return images and targets.
class ClassDataset(Dataset):
    def __init__(self, root, transform=None, demo=False):                
        self.root = root
        self.transform = transform
        self.demo = demo # Use demo=True if you need transformed and original images (for example, for visualization purposes)
        self.imgs_files = sorted(os.listdir(os.path.join(root, "images")))
        self.annotations_files = sorted(os.listdir(os.path.join(root, "annotations")))
    
    def __getitem__(self, idx):
        img_path = os.path.join(self.root, "images", self.imgs_files[idx])
        annotations_path = os.path.join(self.root, "annotations", self.annotations_files[idx])
        name = self.imgs_files[idx]

        img_original = cv2.imread(img_path)
        img_original = cv2.cvtColor(img_original, cv2.COLOR_BGR2RGB)        
        
        with open(annotations_path) as f:
            data = json.load(f)
            bboxes_original = data['bboxes']
            keypoints_original = data['keypoints']
            
            # All objects are glue tubes
            bboxes_labels_original = ['Glue tube' for _ in bboxes_original]            

        if self.transform:   
            # Converting keypoints from [x,y,visibility]-format to [x, y]-format + Flattening nested list of keypoints            
            # For example, if we have the following list of keypoints for three objects (each object has two keypoints):
            # [[obj1_kp1, obj1_kp2], [obj2_kp1, obj2_kp2], [obj3_kp1, obj3_kp2]], where each keypoint is in [x, y]-format            
            # Then we need to convert it to the following list:
            # [obj1_kp1, obj1_kp2, obj2_kp1, obj2_kp2, obj3_kp1, obj3_kp2]
            keypoints_original_flattened = [el[0:2] for kp in keypoints_original for el in kp]
            
            # Apply augmentations
            transformed = self.transform(image=img_original, bboxes=bboxes_original, bboxes_labels=bboxes_labels_original, keypoints=keypoints_original_flattened)
            img = transformed['image']
            bboxes = transformed['bboxes']
            
            # Unflattening list transformed['keypoints']
            # For example, if we have the following list of keypoints for three objects (each object has two keypoints):
            # [obj1_kp1, obj1_kp2, obj2_kp1, obj2_kp2, obj3_kp1, obj3_kp2], where each keypoint is in [x, y]-format
            # Then we need to convert it to the following list:
            # [[obj1_kp1, obj1_kp2], [obj2_kp1, obj2_kp2], [obj3_kp1, obj3_kp2]]
            keypoints_transformed_unflattened = np.reshape(np.array(transformed['keypoints']), (-1,2,2)).tolist()

            # Converting transformed keypoints from [x, y]-format to [x,y,visibility]-format by appending original visibilities to transformed coordinates of keypoints
            keypoints = []
            for o_idx, obj in enumerate(keypoints_transformed_unflattened): # Iterating over objects
                obj_keypoints = []
                for k_idx, kp in enumerate(obj): # Iterating over keypoints in each object
                    # kp - coordinates of keypoint
                    # keypoints_original[o_idx][k_idx][2] - original visibility of keypoint
                    obj_keypoints.append(kp + [keypoints_original[o_idx][k_idx][2]])
                keypoints.append(obj_keypoints)
        
        else:
            img, bboxes, keypoints = img_original, bboxes_original, keypoints_original        
        
        # Convert everything into a torch tensor        
        bboxes = torch.as_tensor(bboxes, dtype=torch.float32)       
        target = {}
        target["boxes"] = bboxes
        target["labels"] = torch.as_tensor([1 for _ in bboxes], dtype=torch.int64) # all objects are glue tubes
        target["image_id"] = torch.tensor([idx])
        target["area"] = (bboxes[:, 3] - bboxes[:, 1]) * (bboxes[:, 2] - bboxes[:, 0])
        target["iscrowd"] = torch.zeros(len(bboxes), dtype=torch.int64)
        target["keypoints"] = torch.as_tensor(keypoints, dtype=torch.float32)        
        img = F.to_tensor(img)
        
        bboxes_original = torch.as_tensor(bboxes_original, dtype=torch.float32)
        target_original = {}
        target_original["boxes"] = bboxes_original
        target_original["labels"] = torch.as_tensor([1 for _ in bboxes_original], dtype=torch.int64) # all objects are glue tubes
        target_original["image_id"] = torch.tensor([idx])
        target_original["area"] = (bboxes_original[:, 3] - bboxes_original[:, 1]) * (bboxes_original[:, 2] - bboxes_original[:, 0])
        target_original["iscrowd"] = torch.zeros(len(bboxes_original), dtype=torch.int64)
        target_original["keypoints"] = torch.as_tensor(keypoints_original, dtype=torch.float32)        
        img_original = F.to_tensor(img_original)

        if self.demo:
            return img, target, img_original, target_original
        else:
            return img, target, name
    
    def __len__(self):
        return len(self.imgs_files)


# Transform the training images and print the output
dataset = ClassDataset(KEYPOINTS_FOLDER_TRAIN, transform=train_transform(), demo=True)
data_loader = DataLoader(dataset, batch_size=5, shuffle=True, collate_fn=collate_fn)

iterator = iter(data_loader)
batch = next(iterator)

# print("Original targets:\n", batch[3], "\n\n")
# print("Transformed targets:\n", batch[1])

# Visualise the transformed images
keypoints_classes_ids2names = {0: 'Hook', 1: 'Tail'}

def visualize(image, bboxes, keypoints, image_original=None, 
    bboxes_original=None, keypoints_original=None, out_file=None):
    fontsize = 18

    for bbox in bboxes:
        start_point = (bbox[0], bbox[1])
        end_point = (bbox[2], bbox[3])
        image = cv2.rectangle(image.copy(), start_point, end_point, (0,255,0), 2)
    
    for kps in keypoints:
        for idx, kp in enumerate(kps):
            image = cv2.circle(image.copy(), tuple(kp), 5, (255,0,0), 10)
            image = cv2.putText(image.copy(), " " + keypoints_classes_ids2names[idx], tuple(kp), cv2.FONT_HERSHEY_SIMPLEX, 1, (255,0,0), 2, cv2.LINE_AA)

    if image_original is None and keypoints_original is None:
        plt.figure(figsize=(40,40))
        plt.imshow(image)
        
        if out_file is None:
            plt.show()
        else:
            plt.savefig(out_file)

    else:
        for bbox in bboxes_original:
            start_point = (bbox[0], bbox[1])
            end_point = (bbox[2], bbox[3])
            image_original = cv2.rectangle(image_original.copy(), start_point, end_point, (0,255,0), 2)
        
        for kps in keypoints_original:
            for idx, kp in enumerate(kps):
                image_original = cv2.circle(image_original, tuple(kp), 5, (255,0,0), 10)
                image_original = cv2.putText(image_original, " " + keypoints_classes_ids2names[idx], tuple(kp), cv2.FONT_HERSHEY_SIMPLEX, 1, (255,0,0), 2, cv2.LINE_AA)

        f, ax = plt.subplots(1, 2, figsize=(40, 20))

        ax[0].imshow(image_original)
        ax[0].set_title('Original image', fontsize=fontsize)

        ax[1].imshow(image)
        ax[1].set_title('Transformed image', fontsize=fontsize)
        if out_file is None:
            plt.show()
        else:
            plt.savefig(out_file)


def show_input(i):        
    image = (batch[0][i].permute(1,2,0).numpy() * 255).astype(np.uint8)
    bboxes = batch[1][i]['boxes'].detach().cpu().numpy().astype(np.int32).tolist()

    keypoints = []
    for kps in batch[1][i]['keypoints'].detach().cpu().numpy().astype(np.int32).tolist():
        keypoints.append([kp[:2] for kp in kps])

    image_original = (batch[2][i].permute(1,2,0).numpy() * 255).astype(np.uint8)
    bboxes_original = batch[3][i]['boxes'].detach().cpu().numpy().astype(np.int32).tolist()

    keypoints_original = []
    for kps in batch[3][i]['keypoints'].detach().cpu().numpy().astype(np.int32).tolist():
        keypoints_original.append([kp[:2] for kp in kps])

    visualize(image, bboxes, keypoints, image_original, bboxes_original, keypoints_original)


# Show example images of transformations
# print("Transforming: ", len(batch), "images for display")
# for i in range(0, len(batch)):
#     print("Showing image", i+1, "of", len(batch))
#     show_input(i)

# Train the Keypoint RCNN model

def get_model(num_keypoints, weights_path=None):
    # anchor_generator = AnchorGenerator()
    anchor_generator = AnchorGenerator(sizes=(48, 64, 128, 256, 512), 
        aspect_ratios=(0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0))

    model = torchvision.models.detection.keypointrcnn_resnet50_fpn(weights=None,
                                                                   weights_backbone=torchvision.models.resnet.ResNet50_Weights.DEFAULT,
                                                                   num_keypoints=num_keypoints,
                                                                   num_classes = 2, # Background is the first class, object is the second class
                                                                   rpn_anchor_generator=anchor_generator)


    if weights_path:
        state_dict = torch.load(weights_path)
        model.load_state_dict(state_dict)        
        
    return model

def train_model(model_path):
	device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')
	# print(device)

	dataset_train = ClassDataset(KEYPOINTS_FOLDER_TRAIN, transform=train_transform(), demo=False)
	dataset_test = ClassDataset(KEYPOINTS_FOLDER_TEST, transform=None, demo=False)

	data_loader_train = DataLoader(dataset_train, batch_size=3, shuffle=True, collate_fn=collate_fn)
	data_loader_test = DataLoader(dataset_test, batch_size=2, shuffle=False, collate_fn=collate_fn)

	model = get_model(num_keypoints = 2)
	model.to(device)

	params = [p for p in model.parameters() if p.requires_grad]
	optimizer = torch.optim.SGD(params, lr=0.001, momentum=0.9, weight_decay=0.0005)
	lr_scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=5, gamma=0.3)
	num_epochs = 10

	for epoch in range(num_epochs):
	    train_one_epoch(model, optimizer, data_loader_train, device, epoch, print_freq=1000)
	    lr_scheduler.step()
	    evaluate(model, data_loader_test, device)
	    
	# Save model weights after training
	torch.save(model.state_dict(), model_path)


# Only run the training if the model does not already exist
model_path = os.path.join(data_path, 'model/weights/'+MODEL_NAME)
if(not os.path.exists(model_path)):
    print("Beginning training")
    train_model(model_path)

print("Model exists")

model = get_model(num_keypoints = 2, weights_path=model_path)

# Make predictions on a subset of images
def make_predictions(in_folder, n_images=None):
    device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')

    batch_size = 5 if n_images is None else n_images

    dataset_test = ClassDataset(in_folder, transform=None, demo=False)
    data_loader_test = DataLoader(dataset_test, batch_size=batch_size, shuffle=True, collate_fn=collate_fn)

    iterator = iter(data_loader_test)
    images, targets, names = next(iterator)
    images = list(image.to(device) for image in images)

    print("Loaded: ", len(images), "images")

    # Run the images against the model
    with torch.no_grad():
        model.to(device)
        model.eval()
        output = model(images)

    # Save the bounding box and keypoint locations
    def save_coordinates(bboxes, keypoints, in_file, out_file):
        with open(out_file, 'a') as f:
            f.write("Image_file\tbbox_x1\tbbox_y1\tbbox_x2\tbbox_y2\tlandmark\tkeypoint_x\tkeypoint_y\n")
            for i in range(len(bboxes)):
                for j in range(len(keypoints[i])):
                    f.write(in_file+"\t"+str(bboxes[i][0])+"\t"+str(bboxes[i][1])+"\t"
                        +str(bboxes[i][2])+"\t"+str(bboxes[i][3])+"\t"
                        +keypoints_classes_ids2names[j]+"\t"+str(keypoints[i][j][0])+"\t"
                        +str(keypoints[i][j][1])+"\n")



    # Visualise predictions
    def show_images(i):
        image = (images[i].permute(1,2,0).detach().cpu().numpy() * 255).astype(np.uint8)
        scores = output[i]['scores'].detach().cpu().numpy()

        high_scores_idxs = np.where(scores > 0.1)[0].tolist() # Indexes of boxes with scores > 0.7
        post_nms_idxs = torchvision.ops.nms(output[i]['boxes'][high_scores_idxs], 
        output[i]['scores'][high_scores_idxs], 0.1).cpu().numpy() # Indexes of boxes left after applying NMS (iou_threshold=0.3)

        # Below, in output[0]['keypoints'][high_scores_idxs][post_nms_idxs] and output[0]['boxes'][high_scores_idxs][post_nms_idxs]
        # Firstly, we choose only those objects, which have score above predefined threshold. This is done with choosing elements with [high_scores_idxs] indexes
        # Secondly, we choose only those objects, which are left after NMS is applied. This is done with choosing elements with [post_nms_idxs] indexes

        keypoints = []
        for kps in output[i]['keypoints'][high_scores_idxs][post_nms_idxs].detach().cpu().numpy():
            keypoints.append([list(map(int, kp[:2])) for kp in kps])

        bboxes = []
        for bbox in output[i]['boxes'][high_scores_idxs][post_nms_idxs].detach().cpu().numpy():
            bboxes.append(list(map(int, bbox.tolist())))
            
        visualize(image, bboxes, keypoints, out_file=os.path.join('./outputs/predication_'+str(i)))
        save_coordinates(bboxes, keypoints, in_file=names[i], out_file=os.path.join('./outputs/predication_'+str(i)+".txt"))

    for i in range(0, len(images)):
        print("Exporting image", i+1, "of", len(images))
        show_images(i)

make_predictions(KEYPOINTS_FOLDER_TEST, n_images=7)
