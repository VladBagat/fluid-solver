import cv2
import glob

files = sorted(glob.glob("animation/frame_*.pgm"))

# load wall mask once
mask = cv2.imread("animation/wall_mask_00000.pgm", cv2.IMREAD_GRAYSCALE)
if mask is not None:
    mask = cv2.resize(mask, (512, 512), interpolation=cv2.INTER_NEAREST)

for f in files:
    smoke = cv2.imread(f, cv2.IMREAD_GRAYSCALE)
    smoke = cv2.resize(smoke, (512, 512), interpolation=cv2.INTER_NEAREST)

    img = cv2.cvtColor(smoke, cv2.COLOR_GRAY2BGR)

    if mask is not None:
        img[mask > 0] = (255, 0, 0)

    cv2.imshow("smoke + walls", img)

    if cv2.waitKey(30) == 27:
        break

cv2.destroyAllWindows()