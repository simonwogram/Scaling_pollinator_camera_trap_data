Bjerge et al., Accurate detection and identification of insects from camera trap images with deep learning, bioRxiv, 2022

Data used to generate Figure 8. The seasonal occurence of individual detections for each of the nine insect species.
Contains the predictions from YOLOv5 best model 1280v5m6 in Table 2.

Name of csv files:
=================
SZ_IP-MonthDate

SZ_IP – Identification of five camera systems: S1_123, S2_146, S3_194, S4_199, S5_187 (Two cameras in each system)
MonthDate – Folder name for where the original image were stored in the system

Format content of cwv files:
===========================
S,C,YYYYMMDD,HHMMSS,Conf,Class,x1,y1,x2,y2,ImageFile

S - System no. 1,2,3,4
C - Camera no. 0,1
YYYYMMDD - Date
HHMMSS - Timestamp
Conf - Class confidence in percentage (0-100)
Class IDs:
1 Coccinellidae septempunctata
2 Apis mellifera
3 Bombus lapidarius
4 Bombus terrestris
5 Eupeodes corolla
6 Episyrphus balteatus
7 Aglais urticae
8 Vespula vulgaris
9 Eristalis tenax

x1,y1,x2,y2 - Bounding box corrdinates of insect found in ImageFile
ImageFile - Path and filename of image
