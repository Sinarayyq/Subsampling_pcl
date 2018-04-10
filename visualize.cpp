#include <visualize.h>

void visualizePointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud1, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud2, std::string label_viewer_window, camera_position camera_pos)
{
	// Window setup
	pcl::visualization::PCLVisualizer viewer(label_viewer_window);
	viewer.setBackgroundColor(0, 0, 0);
	//viewer.addCoordinateSystem(500.0);
	viewer.addCoordinateSystem(10.0, cloud2->at(0).x, cloud2->at(0).y, cloud2->at(0).z, xy);
	viewer.setCameraPosition(0, 0, 1,//double pos_x, double pos_y, double pos_z,                                    
		0, 0, 0,//double view_x, double view_y, double view_z,
		0, 1, 0);//double up_x, double up_y, double up_z, int viewport = 0);

				 // Add the first cloud: size 2, white(default)
	viewer.addPointCloud(cloud1, "cloud1");
	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, "cloud1");
	// Define color "red"
	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> red(cloud2, 255, 0, 0);
	// Add the second cloud: size 3, red (defined)
	viewer.addPointCloud<pcl::PointXYZ>(cloud2, red, "cloud2");
	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "cloud2");

	//Auto recenters the view.
	viewer.resetCamera();

	while (!viewer.wasStopped())
	{
		viewer.spinOnce();
		boost::this_thread::sleep(boost::posix_time::microseconds(100000));
	}
	viewer.close();
}
