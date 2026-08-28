def find_quartets_and_above(min_combinations=4, results_limit=5):
    M = 1  # M represents half the surface area
    results_found = 0
    
    print(f"Searching for the first {results_limit} surface areas with at least {min_combinations} unique cuboids...\n")
    
    # Keep searching until we find the requested number of matching sets
    while results_found < results_limit:
        cuboids = []
        
        a = 1
        while 3 * a * a <= M:
            b = a
            while b * b + 2 * a * b <= M:
                numerator = M - a * b
                denominator = a + b
                
                if numerator % denominator == 0:
                    c = numerator // denominator
                    cuboids.append((a, b, c))
                    
                b += 1
            a += 1
            
        # If this surface area has 4 or more combinations, print it
        if len(cuboids) >= min_combinations:
            results_found += 1
            surface_area = 2 * M
            
            print(f"Result {results_found}: Found {len(cuboids)} cuboids with surface area {surface_area}")
            for i, (a, b, c) in enumerate(cuboids, 1):
                print(f"   Box {i}: {a} x {b} x {c}")
            print("-" * 40)
            
        M += 1
        
    print("Search complete!")

# Run the function to find the first 5 instances of quartets (or larger)
find_quartets_and_above(min_combinations=4, results_limit=5)