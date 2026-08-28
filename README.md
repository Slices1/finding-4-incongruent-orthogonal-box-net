# Finding 4 Incongruent Orthogonal Box Net

the main goal of this project is to find a single 2d polyomino net that can fold into 4 totally unique (incongruent) 3d boxes. as far as i know, a net that folds into 4 boxes has never been found.

we already know it's possible for a net to fold into 2 different boxes. here is an example from the archive of what that looks like:

![2 box net 1](archive/net_that_folds_into_2_unique_boxes_img1.jpg)
![2 box net 2](archive/net_that_folds_into_2_unique_boxes_img2.jpg)

but scaling this up to finding 4 boxes is a massive computational challenge. the project is broken down into three main components:

1. finding surface area quartets
   first we need to find 4 unique boxes that all share the exact same surface area, because they cant share a net if the areas are different. i wrote a python script that mathematically brute forces and filters these quartets.

2. the net generator
   this is where most of the work is right now. it takes a 3d box and generates every single valid 2d net for it. it treats the box surface as a graph and uses a depth first search to build the net square by square. it tracks the 3d to 2d coordinate rotations and backtracks if the paper overlaps itself. we are implementing redelmeier's algorithm to stop it from wasting time checking the exact same net built in a different order.

3. the matcher
   once the generator is fast enough to scale to bigger boxes, this component will take the millions of generated nets from a quartet of boxes and cross reference them to find a single shared shape.

here is the net generator proof of concept working on a standard 1x1x1 cube. it successfully finds all 11 known nets:

![net 1](net_generator_proof_of_concept_images/cube_net_01.png)
![net 6](net_generator_proof_of_concept_images/cube_net_06.png)
![net 11](net_generator_proof_of_concept_images/cube_net_11.png)