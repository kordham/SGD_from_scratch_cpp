/*
Final Assignment, 'The Efficient Setting'
The University of Manchester
School of Physics and Astronomy

Mohammad Kordzanganeh
10137275


This program takes in data from GBplaces.csv and computes
the best number of hubs that would support a network that one
might need to access (send supplies to) all of these Great British places.

Assumptions made:
1. A single package sent serves only one person of the population of a place.
2. All cities are equally interested in our product, same with towns, however towns are more interested,
-because there are no third-party distributors in towns.
3. The delivery vehicles are drones, they move in straight lines and can move only between a place and a hub or two hubs.
4. The drones need to travel between hubs at a median frequency of their trips to towns or cities.

Formula: cos(distance_from_place_to_hub) = cos(latitude_of_hub) cos(latitude_of_place)cos(longitude_of place - longitude of hub)
+sin(latitude of place) sin(latitude of hub)

this can be said to approximately equal ¬ 1- (distance from place to hub)^2 because cos(x)¬1-x^2

Therefore we can find the sum of x^2s and minimise that sum, I have employed this tactic in junction with optimisation for better time efficiency
*/



#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>


#define PI 3.14159265
#define Radius_Earth 6378
//Ratio of the interest of the people living in towns in our product w.r.t that of the cities (regular figure in UK is 1.2 from Online Shopping of towns vs cities)
#define Town_City_Ratio 1.2

using namespace std;


// Structs :: BEGIN


// Position represents a latitude and longitude on the surface of the earth
struct position {

	double latitude=0, longitude=0;
	//addition of two positions
	position addition(position somewhere_else) {
		position result;
		result.longitude = latitude + somewhere_else.latitude;
		result.latitude = longitude + somewhere_else.longitude;
		return result;
	}
	//multiplication of a position by a number, mainly used to add the steps
	position multiplication(double number) {
		position result;
		result.latitude = latitude * number;
		result.longitude = longitude * number;
		return result;
	}
	//setting the position in one line, instead of using two lines
	void set_position(double latitude_temp, double longitude_temp) {
		latitude = latitude_temp;
		longitude = longitude_temp;
	}
	//check if this position is equivalent with another position, returning a boolean, because the '==' operator doesn't work
	bool equals(position somewhere_else) {
		if (somewhere_else.latitude == latitude && somewhere_else.longitude == longitude) {
			return true;
		}
		return false;
	}
};

//a structure storing every line of the GBplaces
struct place_spec {

	string name;
	//the following stores true if it's a city and false otherwise
	bool city;
	double population;
	position position;
	//weight() gives the weight of the place for the purposes of maximising the profit
	double weight() {
		double weight;
		//favouring more populated places
		weight = population;
		//favouring towns
		if (!city) {
			weight *= Town_City_Ratio;
		}
		return weight;
	}
};

//best output stores the outcome of an iteration
struct best_output {
	position coords;
	double distance;
};

//grad_hat_constants stores the constants associated with the gradient of the approximated function
struct grad_hat_constants {
	double cos_fi_cos_l = 0, sin_l_cos_fi = 0, sin_fi = 0;
};

// a result, which includes the total distance of all the hubs to their corresponding places as well as each-other +
// the positions of the hubs themselves
struct result {
	double total_distance=0;
	vector<position> hubs;
};

//Structs :: END


//implementing the bubble sort algorithm to sort the file w.r.t latitude
vector<place_spec> bubble_sort(vector<place_spec> data) {

	for (int i = 0; i < data.size() - 1; i++) {
		// Last i elements are already in place    
		for (int j = 0; j < data.size() - i - 1; j++) {
			//check if jth latitude is greater than j+1, if so swap them around and keep
			//doing this until this doesn't happen again
			if (data[j].position.latitude > data[j + 1].position.latitude) {
				place_spec temp = data[j];
				data[j] = data[j + 1];
				data[j + 1] = temp;
			}
		}
	}
	return data;
}

//fetching the data
vector<place_spec> fetch_data(string file_read_name) {

	//declare 2D vector matrix, a string that will contain the lines and a boolean value
	vector <place_spec> data;
	place_spec GB_place;
	string line;
	//this boolean value indicates whether an error message has already been shown once
	bool inconsistent_msg_seen = false;

	//insert the file
	ifstream file_pointer(file_read_name);

	//Was the insertion sucessful?
	if (file_pointer.is_open()) {

		//declare the number of items to be the maxmimum possible value so as to be referenced later
		int number_items = -1;
		int comma, line_number = 1;

		//while-loop below iterates until the end of file is reached
		while (!file_pointer.eof()) {
			//fetch each line
			getline(file_pointer, line);

			//don't read the first line
			if (line_number != 1) {

				//count function returns the number of occurences of a character, ',' in this case, in a string.
				//note that there is one more item than the number of commas, as they split the items.
				number_items = static_cast<int>(count(line.begin(), line.end(), ',') + 1);

				// need a temporary vector to hold push_back our data from, comma starts from -1, because 0 would be an actual possible value
				vector <double> temp_matrix;
				comma = -1;

				//for a given line, iterate through every item and store it in the temporary vector, temp_matrix
				for (int i = 0; i < number_items; i++) {
					// store comma in another variable as it needs to be changed itself but its current value will be needed too
					int temp_comma = comma;
					//find where the comma is and call it comma
					comma = line.find(',', comma + 1);
					//choose the substring starting from right after previous comma and ending right before the current comma
					string item = line.substr(temp_comma + 1, comma - temp_comma - 1);

					//entering the data from each column to the corresponding place in a place_spec
					switch (i) {
					
					case 0: GB_place.name = item;
					case 1:
						if (item == "City") { GB_place.city = true; }
						else {
							GB_place.city = false;
						}; break;
					case 2: GB_place.population = stod(item); break;
					case 3: GB_place.position.latitude = stod(item) / 180 * PI; break;
					case 4: GB_place.position.longitude = stod(item) / 180 * PI; break;
					}

				}

				//inset the elements on this line into the actual matrix
				data.push_back(GB_place);


			}

			//add the number of line by one to keep track of which row we're visiting
			line_number++;
		}
		//close the file
		file_pointer.close();

		return data;
	}

	// if the opening of the file was unsuccessful, alert the user and exit the program
	else {
		file_pointer.close();
		cout << "Couldn't open and read the file, ' " << file_read_name << " '. Please check the directory again." << endl;
		system("pause");
		exit(0);
	}



}
//fetch data end


//grad_hat find the direction in which the function is increasing the most, so we can then decide to go against that
position grad_hat(position coords,grad_hat_constants constants) {
	position grad_hat;

	double d_by_dlongitude = 0, d_by_dlatitude = 0,grad_modulus=0;
	
	d_by_dlongitude = sin(coords.latitude)*(sin(coords.longitude)*constants.cos_fi_cos_l - cos(coords.longitude)*constants.sin_l_cos_fi);
	d_by_dlatitude = sin(coords.latitude)*(cos(coords.longitude)*constants.cos_fi_cos_l + sin(coords.longitude)*constants.sin_l_cos_fi)
		- cos(coords.latitude)*constants.sin_fi;
	// to find the hat, we need to divide through by the modulus
	grad_modulus = sqrt(pow(d_by_dlatitude, 2) + pow(d_by_dlongitude, 2));
	//set the found positions and divide by the modulus
	grad_hat.set_position(d_by_dlongitude, d_by_dlatitude);
	grad_hat.multiplication(1 / grad_modulus);

	return grad_hat;
}
// grad hat end

//grad hat constants finds the constants needed to calculate the grad_hat
grad_hat_constants grad_summer(vector <place_spec> data) {
	grad_hat_constants summer_output;
	//summing over each value required for the grad value;
	for (size_t i = 0; i < data.size(); i++) {
		summer_output.cos_fi_cos_l += cos(data[i].position.latitude)*cos(data[i].position.longitude)*pow(data[i].weight(), 2);
		summer_output.sin_fi += sin(data[i].position.latitude)*pow(data[i].population,2);
		summer_output.sin_l_cos_fi = sin(data[i].position.longitude)*cos(data[i].position.latitude)*pow(data[i].weight(), 2);
	}
	return summer_output;
}
//end grad_hat_constants


//find the distance between two places
double solo_distance_finder(position place1, position place2) {

	double distance;
	// use the great circle distance formula
	distance = (acos(sin(place1.latitude)*sin(place2.latitude)
		+ cos(place2.latitude)*cos(place1.latitude)*cos(place2.longitude - place1.longitude)))*Radius_Earth;
	return distance;
}
// end solo_distance_finder

//finds the total distance by summing over all places, this can be given a boolean to indicate whether we need the distances weighted or not
double total_distance_finder(vector<place_spec> data, position place,bool weight) {
	double total_distance = 0;
	for (size_t i = 0; i < data.size(); i++) {
		//if weighted each solo distance is multiplied by a weight
		if (weight) {
			total_distance += solo_distance_finder(data[i].position, place)*data[i].weight();
		}
		else {
			total_distance += solo_distance_finder(data[i].position, place);
		}
	}

	return total_distance;
	
}
//total_distance_finder end

//initial_position finds a very good approximation of the initial position using the approximating function
position initial_position(vector<place_spec> data) {
	//declaring numerator and denominator
	double top = 0, bottom = 0;
	position coords;
	
	for (size_t i = 0; i < data.size(); i++) {
		//weighting by squared weight because this approximation assumes cosx ¬ 1-x^2
		bottom += cos(data[i].position.latitude)*cos(data[i].position.longitude)*pow(data[i].weight(), 2);
		top += cos(data[i].position.latitude)*sin(data[i].position.longitude)*pow(data[i].weight(), 2);
	}
	//finding the longitude
	coords.longitude = atan2(top, bottom);
	//nullifying the top and bottom for future use
	top = 0; bottom = 0;

	for (size_t i = 0; i < data.size(); i++) {
		
		bottom += cos(data[i].position.latitude)*cos(data[i].position.longitude - coords.longitude)*pow(data[i].weight(), 2);
		top += sin(data[i].position.latitude)*pow(data[i].weight(), 2);


	}
	//finding the latitude
	coords.latitude = atan2(top, bottom);

	return coords;
}
//initial_position end


//the recursor function moves in the negative direction of the gradient hat of the approximated function to get close to the 
best_output recursor(vector<place_spec> data, position coords, double step, grad_hat_constants constants) {
	best_output coords_distance;
	//defining the direction of grad
	position grad_direction = grad_hat(coords, constants)
		, alternative;
	//adding the grad *step to the position, notice it's timsed by -1 to move in the opposite direction of the increase
	alternative = alternative.addition(coords.addition(grad_direction.multiplication(-1 * step)));

	double distance = total_distance_finder(data, coords, true);

	//if the total distance is greater than the distance of the next step, redo the function for the next step
	if (distance > total_distance_finder(data, alternative, true)) {
		coords = alternative;
		coords_distance = recursor(data, coords, step, constants);
	}
	else {
		//else, reduce your step size, to get more precise
		if (step > 0.01*PI/180) {
			coords_distance = recursor(data, coords, step*0.8, constants);
		}//when precise enough, set the coordinates to the returning value
		else {
			coords_distance.coords = coords;
			coords_distance.distance = distance;
		}

	}
    // return that value
	return coords_distance;
}
//recursor end



//iterate is the next step after recursor, becuase recursor is an approximation, this is used to increase accuracy
best_output iterate(vector <place_spec> data, position location, double step) {
	
	best_output coords;
	// find the total distance weighted
	double distance = total_distance_finder(data, location,true);
	position temp, store;
	//we assume a smaller value is not found
	bool smaller_value_found = false;
	// we go through the 8 blocks adjacent and check for a smaller value
	for (int i = -1; i <= 1; i++) {
		for (int j = -1; j <= 1; j++) {
			//skip the same block 
			if (i == 0 && j == 0) {}
			else {
				//if there is a smaller distance somewhere, set the boolean to true and store its points
				temp.longitude = location.longitude + i * step;
				temp.latitude = location.latitude + j * step;
				if (distance > total_distance_finder(data, temp,true)) {
					store = temp;
					smaller_value_found = true;

				}
			}
		}
	}
	//if smaller value was found, reiterate the iterate 
	if (smaller_value_found) {
		coords = iterate(data, store, step); }
	else {
	// if no smaller values found, set the variables and return the best_output
			coords.coords = location; coords.distance = distance;
	}
	return coords;

}
//iterate end



//optimisation, this function takes in the data and tries to use the past functions to produce a best position for a single hub
position optimisation(vector<place_spec> data) {
	
	//declare a place_holder
	best_output resultant;
	//calculate the constants of grad
	grad_hat_constants constants = grad_summer(data);
	bool flag = true;
	// set the minimum disatnce to the maximum possible value
	double min_distance = DBL_MAX;

	//start the recursor from the initial positions and use 0.1 degree intervals
	resultant = recursor(data, initial_position(data), 0.1 / 180 * PI, constants);
	//then smoothen the sharp edges using iterate with 0.01 precision
	resultant = iterate(data,resultant.coords, 0.01/180*PI);

	//return the coordinates
	return resultant.coords;

}
//end optimisation


//n_hubs_solution, repeats the same procedure but for n hubs(takes in the number of hubs)
result n_hubs_solution(vector<place_spec> data,int number_hubs) {
	vector<vector<place_spec>> hubs;
	vector <position> hub_positions;
	vector <place_spec> filler;
	position filler_position;
	//need to initialise these vectors first
	for (int i = 0; i < number_hubs; i++) {
		hubs.push_back(filler);
		hub_positions.push_back(filler_position);
	}
	//this loop assigns different places to different hubs, by assigning them in order of latitude, as the parse data 
	//argument is sorted in that order
	int counter = 0;
	for (int i = 0; i < data.size(); i += data.size() / number_hubs) {
		
		for (int j = 0; j < data.size() / number_hubs-1; j++) {
			//give all the remaining places to the last hub
			if ((i+j)>=(int)((data.size()/number_hubs))*number_hubs) {
				counter--;
			}
			if ((i + j )< data.size()) {
				
				hubs[counter].push_back(data[i + j]);
			}
		}
		
		counter++;

	}

	//clustering the data, the following sequentially changes the assignment of the data to the hubs
	int count = 0;
	position temp;
	bool ended_cycle = false;
	//iterates 10 times
	while (count < 10) {
		count++;
		
		for (int k = 0; k < number_hubs; k++) {
			//find the position of the hubs
			hub_positions[k] = optimisation(hubs[k]);
			//if ALL hubs are remaining the same break out of the loop
			if (hub_positions[k].equals(temp)) {
				ended_cycle = true;
			}
			else {
				ended_cycle = false;
			}
			//remove the contents of the hubs
			hubs[k].clear();
		}
		//break if no change
		if (ended_cycle) {
			break;
		}

		int index;

		//assign every point to its nearest hub
		for (size_t i = 0; i < data.size(); i++) {
			double min_distance = DBL_MAX;
			for (int k = 0; k < number_hubs; k++) {
				double temp_distance = solo_distance_finder(data[i].position, hub_positions[k]);
				//find the nearest hub and assign its index 'k' to index and distance to min_distance
				if (temp_distance < min_distance) {
					min_distance = temp_distance;
					index = k;

				}
			}
			//add that place to the index hub
			hubs[index].push_back(data[i]);

		}
		count++;
		//repeat this to slowly get better results
	}

	result all_hubs;
	//finding the distance of all hubs to their corresponding places without weighting
	for (int i = 0; i < number_hubs; i++) {
		all_hubs.hubs.push_back(hub_positions[i]);
		all_hubs.total_distance += total_distance_finder(hubs[i], hub_positions[i], false);
	}
	// The handshake problem, every hub goes to all the other hubs, which means the path between any two hubs is taken twice
	double twice_distance_hubs = 0;
	for (int i = 0; i < number_hubs; i++) {
		for (int j = 0; j < i; j++) {
			twice_distance_hubs += solo_distance_finder(hub_positions[i], hub_positions[j]);
		}
	}
	//halving the distance because of the mentioned problem
	all_hubs.total_distance += twice_distance_hubs / 2;
	//returning it
	return all_hubs;

}
//n_hubs_solution end


//MAIN :: BEGIN
int main() {
	

	//we fetch and sort the data from the file GBplaces.csv
	vector<place_spec> data = bubble_sort(fetch_data("GBplaces.csv"));
	

	double min_distance = DBL_MAX,temp_distance;
	vector<position> temp_hub_holder;


	//number of hubs is n
	int number_of_hubs = 20,
		index;
	//this for loop checks every possibility of 1, 2,..., n hubs to see which one is the most efficient
	for (int i = 1; i <= number_of_hubs; i++) {
		
		result temp_result = n_hubs_solution(data, i);
		temp_distance = temp_result.total_distance;
		//find which one poses the least distance in total
		if (temp_distance < min_distance) {
			//minimum distance uses the same algorithm as in optimiser
			min_distance = temp_distance;
			index = i;
			temp_hub_holder = temp_result.hubs;
		}
		//all the distance for each of the settings
		cout << "The cummulative distance for a setting with " << i << " hub(s) is = " << temp_distance << endl;
	}
	//most efficient?
	cout << "Therefore the most efficient way is when number of hubs : " << index << endl << endl;

	//where each hub is situated
	for (int i = 0; i < index; i++) {
		//print each of them in degrees
		cout << "Hub Number (" << i+1 << ") is situated at latitude "
			<< temp_hub_holder[i].latitude / PI * 180 << " and longitude " << temp_hub_holder[i].longitude / PI * 180 <<endl;
	}
	

	system("pause");

	return 0;
}
//MAIN :: END