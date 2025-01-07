# RLML: Reflector-based LiDAR Mapping and Localization

---

This project uses **2D LiDAR and odometry**, with **reflector assistance**, to achieve **mapping and localization**.  
Key features of the project include:
- **Mapping** based on LiDAR scans and reflectors
- **Relocalization** using the generated map
- **Real-time localization**, with **odometry smoothing** to provide a stable trajectory for more reliable control

---
## Demo data
Download the demo data in the following link:
https://drive.google.com/drive/folders/1j5oe1IJGGBV9VHBFqLp63lBu8AV0TvW1?usp=sharing

## How to Build and Run

### 1. Install Dependencies
Before starting, please ensure the following dependencies are installed:
- **ROS1**
- **PCL**
- **Cholmod**

### 2. Build g2o
```bash
cd Thirdparty/g2o
mkdir build && cd build
cmake ..
make
```

### 3. Build the Main Project


```bash
cd ../../..
mkdir build && cd build
cmake ..
make
```
### 4. Run Mapping


```bash
export ROS_PACKAGE_PATH=<your project root>:$ROS_PACKAGE_PATH
roslaunch RLML mapping.launch
```
Download the **demo dataset**  from the provided link and run:

```bash
rosbag play --clock mapping.bag
```
Once the data finishes playing, type `s`  and press Enter in the pop out window of `mapping.launch` to save the map.
### 5. Run Localization


```bash
export ROS_PACKAGE_PATH=<your project root>:$ROS_PACKAGE_PATH
roslaunch RLML localization.launch
```
Download the **demo dataset**  from the provided link and run:

```bash
rosbag play --clock localization.bag
```
After the data finishes playing, type `q`  and press Enter in the pop out window of  `localization.launch` to exit the program.

---

## TODO List

- [ ] Add support for ROS2
- [ ] Add a docker support


## Acknowledgements

This project was inspired by the following open-source projects:

- [g2o](https://github.com/RainerKuemmerle/g2o): A General Framework for Graph Optimization.
- [ares_slam](https://github.com/ningwang1028/ares_slam): A LiDAR-based SLAM system.

## License
This project is licensed under the **MIT License** .
For more details, see the [LICENSE]()  file.

---

