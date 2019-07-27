# NFL Big Data Bowl 2018

Contact: Mitchell Kinney, kinne174@umn.edu

My entry into the NFL Big Data Bowl and subsequent paper currently under review with 
the Journal of Quantitative Analysis in Sports. My goal was to classify
routes and use those classifications to predict if certain route 
combinations proved more successful than others when running a
passing play in the NFL. 

I used a nearest neighbor metric in Euclidean distance to match routes
run by receivers to a manually defined route tree. Then I gathered 
spatial information about the players while the pass play was being
run such as the separation the intended receiver had when the ball was
thrown. Finally based on my own definition of success, I used the 
classified routes to predict the spatial information as a stepping
stone to then predict success of the play. 

The models used were a Gibbs sampling set up to predict the spatial
information and then a penalized logistic regression to predict
success. 

An example of my route classification for a random game:

![route classification](https://github.com/kinne174/BigDataBowl_2018/blob/master/pictures/all_routes.PNG)

