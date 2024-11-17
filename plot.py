import matplotlib.pyplot as plt

# Read data from file
x_coords = []
y_coords = []

with open("orbit_data.txt", "r") as file:
    for line in file:
        x, y = map(float, line.split())
        x_coords.append(x)
        y_coords.append(y)

# Plot the orbit
plt.figure(figsize=(8, 8))
plt.plot(x_coords, y_coords, label="Orbit Path")
# plt.scatter(0, 0, color="red", label="Central Body")  # Mark the central body
plt.xlabel("x position (meters)")
plt.ylabel("y position (meters)")
plt.title("Elliptical Orbit Plot")
# plt.legend()
plt.grid(True)
plt.axis("equal")
plt.show()
