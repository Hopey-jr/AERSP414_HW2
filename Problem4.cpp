
// AERSP 424, Homework 2 Problem 4
/*This code plots the orbital energy equation, specifically energy vs velocity over time using openCV. The first plot is 
a normal plot, but the second plot is the data from the CSV file that was written and read in this code*/ 

/*Sources 
https://stackoverflow.com/questions/15955305/find-maximum-value-of-a-cvmat
https://www.geeksforgeeks.org/python-opencv-cv2-puttext-method/
https://www.geeksforgeeks.org/how-to-read-data-from-csv-file-to-a-2d-array-in-cpp/
https://www.geeksforgeeks.org/csv-file-management-using-c/

*/

#include <iostream>
#include <opencv2/opencv.hpp>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cmath>


int main()
{
    // create all of the variables
    int width = 800;     // plot width
    int height = 600;   // plot height
    double mu = 398600e9; // gravitation constant in m^3/s^2
    double r0 = 8000e3; //  initial distance from spacecraft to center of earth in meters
    double v0 = sqrt(mu/r0); // initial velocity in m/s
    double dt = 0.1; // time step
    double duration = 5000; // duration of simulation in seconds

    // create widget for the plot
    cv::Mat image = cv::Mat::zeros(height, width, CV_8UC3);

    // set initial position and velocity to be r0, v0
    double r = r0;
    double v = v0;

    // create vectors to hold data gathered from simulation
    std::vector<double> time;
    std::vector<double> velocity;
    std::vector<double> energy;

    // Write data to CSV file - this is for part 2. If the csv file has not been created, write it, and 
    // if the file is alredy in existence, overwrite the file with new content from simulation
    std::ofstream file("data.csv");
    if (!file.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        return -1;
    }
    // Write header
    file << "Time,Velocity,Energy\n";  

    // generate points to plot - from 0 to end of duration incremented by dt
    for (double t = 0; t < duration; t += dt){
        
        // Energy equation 
         double orb_energy = 0.5*v*v - mu/r;
       
       // store data in vectors
        time.push_back(t);
        velocity.push_back(v);
        energy.push_back(orb_energy);

        // Simple numerical integration (Euler method for velocity and position)
        double acceleration = mu / (r * r); 

        // Update velocity and position
        v += acceleration * dt;  
        r += v * dt;             

        // Prevent r from becoming less than the radius of the earth
        if (r < 6378e3) {
            r = 6378e3;}

        // Write data to CSV file with correct headers
        file << t << "," << v << "," << orb_energy << "\n";
     }
        
    // close the csv file
    file.close();

    // find max and min values in velocity and energy vectors -  use later for scaling the plot
    double max_velocity = *std::max_element(velocity.begin(), velocity.end());
    double min_velocity = *std::min_element(velocity.begin(), velocity.end());
    double max_energy = *std::max_element(energy.begin(), energy.end());
    double min_energy = *std::min_element(energy.begin(), energy.end());

    // Calculate scaling factors based on min and max values for velocity and energy
    double scale_v = (width - 100) / (max_velocity - min_velocity);  // scale for x (velocity)
    double scale_e = (height - 100) / (max_energy - min_energy);      // scale for y (energy)

    // Create the plot
    for (size_t i = 1; i < velocity.size(); ++i ){
 
         // Calculate x and y coordinates for the two points to create the line
         // this section scales the v and e values and shifts by 50
         int v1 = static_cast<int>((velocity[i - 1] - min_velocity) * scale_v) + 50; 
         int v2 = static_cast<int>((velocity[i] - min_velocity) * scale_v) + 50;      

         int E1 = static_cast<int>((max_energy - energy[i - 1]) * scale_e) + 50;    
         int E2 = static_cast<int>((max_energy - energy[i]) * scale_e) + 50;          

         // put data on plot by creating line at start and end points
         cv::line(image, cv::Point(v1, E1), cv::Point(v2, E2), cv::Scalar(0, 255, 0), 2);
    }
    // PLOT DETAILS
    // x axis and label
    cv::line(image, cv::Point(50, height - 50), cv::Point(width - 50, height - 50), cv::Scalar(255, 255, 255), 2); 
    cv::putText(image, "Velocity (m/s)", cv::Point(width-480, height - 10), cv::FONT_HERSHEY_SIMPLEX, 0.7, cv::Scalar(255, 255, 255), 2);

    // y axis and label
    cv::line(image, cv::Point(50, height - 50), cv::Point(50, 50), cv::Scalar(255, 255, 255), 2);  
    cv::putText(image, "Energy (J)", cv::Point(10, 30), cv::FONT_HERSHEY_SIMPLEX, 0.7, cv::Scalar(255, 255, 255), 2);
    
    // Add gridlines and ticks
    int tick_size = 10;  
    int grid_spacing = 50;  

    // Vertical grid lines
    for (int i = 50; i < width - 50; i += grid_spacing) {
        cv::line(image, cv::Point(i, 50), cv::Point(i, height - 50), cv::Scalar(255, 255, 255), 1); 
    }
    // Horizontal grid lines
    for (int i = height - 50; i > 50; i -= grid_spacing) {
        cv::line(image, cv::Point(50, i), cv::Point(width - 50, i), cv::Scalar(255, 255, 255), 1); 
    }

    // X-axis tick marks
    for (int i = 50; i < width - 50; i += grid_spacing) {
        cv::line(image, cv::Point(i, height - 50), cv::Point(i, height - 50 + tick_size), cv::Scalar(255, 255, 255), 1);

        // Calculate the velocity value corresponding to the X position - for tick marks
        double velocity_value = min_velocity + (i - 50) * (max_velocity - min_velocity) / (width - 100);
        // convert to integer to handle plot 
        std::string velocity_label = std::to_string(static_cast<int>(velocity_value));  
        cv::putText(image, velocity_label, cv::Point(i - 15, height - 40), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(255, 255, 255), 1);
    }

    // Y-axis tick marks
    for (int i = 50; i < height - 50; i += grid_spacing) {
        cv::line(image, cv::Point(50 - tick_size, i), cv::Point(50, i), cv::Scalar(255, 255, 255), 1);

        // Calculate the energy value corresponding to the y position - for tick marks
        double energy_value = min_energy + (height - 50 - i) * (max_energy - min_energy) / (height - 100);
        // convert to integer to handle plot 
        std::string energy_label = std::to_string(static_cast<int>(energy_value));  // Integer value for readability
        cv::putText(image, energy_label, cv::Point(10, i + 5), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(255, 255, 255), 1);
    }

    // add title and show plot
    cv::putText(image, "Energy vs Velocity", cv::Point(140, 30), cv::FONT_HERSHEY_SIMPLEX, 1, cv::Scalar(255, 255, 255), 2);
    cv::imshow("Energy vs Velocity", image);
    cv::waitKey(0);

    // print out data for debugging
     //for (size_t i = 1; i < velocity.size(); ++i) {
     //std::cout << "Time: " << time[i] << "s, Velocity: " << velocity[i] << " m/s, Energy: " << energy[i] << " J" << std::endl;
    //}

// PART 2 OF HW - read data from CSV and plot

    // Open CSV file in read mode
    std::ifstream inputfile("data.csv");
    if (!inputfile.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        return -1;
    }

    // Store data in CSV as time, velocity, energy
    std::vector<double> Time_csv;
    std::vector<double> Velocity_csv;
    std::vector<double> Energy_csv;

    // Read the header line ("Time, Velocity, Energy")
    std::string line;
    std::getline(inputfile, line);  
    
    // Read the rest of the lines containing data
    while (std::getline(inputfile, line)) {
        std::stringstream ss(line);
        // create strings for each 
        std::string time_str, v_str, energy_str;

        // read time
        std::getline(ss, time_str, ',');  
        // read velocity    
        std::getline(ss, v_str, ','); 
        // read energy      
        std::getline(ss, energy_str, ','); 

        // Convert the string data to double and store in vectors to use for plotting
        Time_csv.push_back(std::stod(time_str));
        Velocity_csv.push_back(std::stod(v_str));
        Energy_csv.push_back(std::stod(energy_str));
    }

    // Close the file after reading
    file.close();  

    // Plotting the CSV data

    // initialize new width and height
    int width2 = 800; 
    int height2 = 600;

    // create widget for the plot
    cv::Mat plot = cv::Mat::zeros(height2, width2, CV_8UC3);

    // Find the min and max values of velocity and energy vectors - scaling the plot
    double min_v = *std::min_element(Velocity_csv.begin(), Velocity_csv.end());
    double max_v = *std::max_element(Velocity_csv.begin(), Velocity_csv.end());
    double min_E = *std::min_element(Energy_csv.begin(), Energy_csv.end());
    double max_E = *std::max_element(Energy_csv.begin(), Energy_csv.end());

    // Scaling factors - subtract 100 and normalize for correct scale
    double scale_x = (width2 - 100) / (max_v - min_v); 
    double scale_y = (height2 - 100) / (max_E - min_E);

    // Create axis
    // x axis
    cv::line(plot, cv::Point(50, height2 - 50), cv::Point(width2 - 50, height2 - 50), cv::Scalar(255, 255, 255), 2);  
    // y axis
    cv::line(plot, cv::Point(50, height2 - 50), cv::Point(50, 50), cv::Scalar(255, 255, 255), 2);  

    // Plot data 
    for (size_t i = 1; i < Time_csv.size(); ++i) {

        // convert back to integer to handle plotting. Scale by subtracting min value 
        // and then multiplying by scale. Then add 50 to shift
        int x1 = static_cast<int>((Velocity_csv[i - 1] - min_v) * scale_x) + 50;
        int x2 = static_cast<int>((Velocity_csv[i] - min_v) * scale_x) + 50;

        int y1 = static_cast<int>((max_E - Energy_csv[i - 1]) * scale_y) + 50;
        int y2 = static_cast<int>((max_E - Energy_csv[i]) * scale_y) + 50;

        // put data on plot by creating line
        cv::line(plot, cv::Point(x1, y1), cv::Point(x2, y2), cv::Scalar(0, 255, 0), 2);
    }

    // PLOT DETAILS FOR CSV PLOT
    // initialize new grid and ticks 
    int grid_space = 50;
    int ticks = 10;

    // Add gridlines (vertical and horizontal lines)
    // Vertical grid lines
    for (int i = 50; i < width2 - 50; i += grid_space) {
        cv::line(plot, cv::Point(i, 50), cv::Point(i, height2- 50), cv::Scalar(255, 255, 255), 1); 
    }
    // Horizontal grid lines
    for (int i = height2 - 50; i > 50; i -= grid_space) {
        cv::line(plot, cv::Point(50, i), cv::Point(width2 - 50, i), cv::Scalar(255, 255, 255), 1); 
    }

    // X axis tick marks and scale
    for (int i = 50; i < width2 - 50; i += grid_space) {
        cv::line(plot, cv::Point(i, height2 - 50), cv::Point(i, height2 - 50 + ticks), cv::Scalar(255, 255, 255), 1);

        // Calculate the velocity value corresponding to the X position
        double velocity_new = min_v + (i - 50) * (max_v - min_v) / (width2 - 100);
        // convert to integer to handle plot
        std::string velocity_label2 = std::to_string(static_cast<int>(velocity_new));  
        cv::putText(plot, velocity_label2, cv::Point(i - 15, height2 - 40), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(255, 255, 255), 1);
    }

    // Y axis tick marks and scale
    for (int i = 50; i < height2 - 50; i += grid_space) {
        cv::line(plot, cv::Point(50 - ticks, i), cv::Point(50, i), cv::Scalar(255, 255, 255), 1);
        
        // Calculate the energy value corresponding to the Y position      
        double energy_new = min_E + (height2 - 50 - i) * (max_E - min_E) / (height2 - 100);
        // convert to integer to handle plot
        std::string energy_label2 = std::to_string(static_cast<int>(energy_new));  
        cv::putText(plot, energy_label2, cv::Point(10, i + 5), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(255, 255, 255), 1);
    }

    // Add labels to axes
    cv::putText(plot, "Velocity (m/s)", cv::Point(width2 - 480, height2 - 10), cv::FONT_HERSHEY_SIMPLEX, 0.7, cv::Scalar(255, 255, 255), 2);
    cv::putText(plot, "Energy (J)", cv::Point(10, 30), cv::FONT_HERSHEY_SIMPLEX, 0.7, cv::Scalar(255, 255, 255), 2);
    cv::putText(plot, "CSV 2D Plot of Energy vs Velocity", cv::Point(140, 30), cv::FONT_HERSHEY_SIMPLEX, 1, cv::Scalar(255, 255, 255), 2);

    // display image
    cv::imshow("CSV Data Plot", plot);
    cv::waitKey(0); 

return 0;
}
