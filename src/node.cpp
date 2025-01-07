#include "RLML_ros.h"


int main(int argc, char** argv)
{
    ros::init(argc, argv, "RLML_node");
    ros::console::set_logger_level(ROSCONSOLE_DEFAULT_NAME, ros::console::levels::Info);


    RLMLRos mapper;

    ros::spin();
    return 0;
}
