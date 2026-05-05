import cv2
import glob

files = sorted(glob.glob("animation/frame_*.pgm"))

for f in files:
    smoke = cv2.imread(f, cv2.IMREAD_GRAYSCALE)
    smoke = cv2.resize(smoke, (512, 512), interpolation=cv2.INTER_NEAREST)

    img = cv2.cvtColor(smoke, cv2.COLOR_GRAY2BGR)

    cv2.imshow("smoke", img)

    if cv2.waitKey(30) == 27:
        break

cv2.destroyAllWindows()
