#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>

typedef struct
{
    double x;
    double y;
} Point;

double model(double x, double a, double b, double c)
{
    return a * sin(b * x + c);
}

double mean_squared_error(Point points[], int n, double a, double b, double c)
{
    double error = 0.0;
    for (int i = 0; i < n; i++)
    {
        double diff = model(points[i].x, a, b, c) - points[i].y;
        error += diff * diff;
    }
    return error / n;
}

void compute_gradients(Point points[], int n, double a, double b, double c, double *grad_a, double *grad_b, double *grad_c)
{
    *grad_a = 0.0;
    *grad_b = 0.0;
    *grad_c = 0.0;

    for (int i = 0; i < n; i++)
    {
        double x = points[i].x;
        double error = model(x, a, b, c) - points[i].y;

        *grad_a += 2 * error * sin(b * x + c) / n;         // Parcialni derivace vzhledem k a
        *grad_b += 2 * error * a * x * cos(b * x + c) / n; // Parcialni derivace vzhledem k b
        *grad_c += 2 * error * a * cos(b * x + c) / n;     // Parcialni derivace vzhledem k c
    }
}

void gradient_descent_sine_wave(Point points[], int n, double *a, double *b, double *c, double rate, int iterations)
{

    for (int i = 0; i < iterations; i++)
    {
        double grad_a, grad_b, grad_c;
        compute_gradients(points, n, *a, *b, *c, &grad_a, &grad_b, &grad_c);

        *a -= rate * grad_a;
        *b -= rate * grad_b;
        *c -= rate * grad_c;

        if (i % 100 == 0)
        {
            double error = mean_squared_error(points, n, *a, *b, *c);
            printf("Iterace %d: a = %f, b = %f, c = %f, Chyba: %f\n", i, *a, *b, *c, error);
        }
    }
}

void iterate_sine_wave(Point points[], int n, double *best_a, double *best_b, double *best_c,
                       double a_min, double a_max, double a_step,
                       double b_min, double b_max, double b_step,
                       double c_min, double c_max, double c_step)
{
    double min_error = DBL_MAX;

    for (double a = a_min; a <= a_max; a += a_step)
    {
        for (double b = b_min; b <= b_max; b += b_step)
        {
            for (double c = c_min; c <= c_max; c += c_step)
            {
                double error = mean_squared_error(points, n, a, b, c);

                if (error < min_error)
                {
                    min_error = error;
                    *best_a = a;
                    *best_b = b;
                    *best_c = c;
                }
            }
        }
    }

    printf("Chyba: %f\n", min_error);
}

int main()
{
    // To use, rename the list to points
    // Io
    Point Io[] = {
        {0, 1.88},
        {0.1944, 0.00},
        {0.9549, -1.30},
        {1.9514, 0.15},
        {2.9410, 0.99}};
    // Europa
    Point Europa[] = {
        {0, -3.41},
        {0.9549, 3.79},
        {1.9514, 2.33},
        {2.9410, -4.51}};
    // Ganymedes
    Point Ganymedes[] = {
        {0, -5.8},
        {0.9549, -0.43},
        {1.9514, 5.77},
        {2.9410, 7.54}};
    // Callisto
    Point points[] = {
        {0, 13.30},
        {0.9549, 13.50},
        {1.9514, 11.80},
        {2.9410, 8.49},
        {4.7535, 0.00},
        {8.9827, -13.40}};

    int n = sizeof(points) / sizeof(points[0]);

    double rate = 0.001;
    int iterations = 1000;

    double a_min = 0, a_max = 20, a_steps = 1000, a_step = (a_max - a_min) / a_steps;
    double b_min = 0, b_max = 5, b_steps = 1000, b_step = (b_max - b_min) / b_steps;
    double c_min = 0, c_max = 2 * M_PI, c_steps = 1000, c_step = (c_max - c_min) / c_steps;

    double best_a = 0.0, best_b = 0.0, best_c = 0.0;

    iterate_sine_wave(points, n, &best_a, &best_b, &best_c, a_min, a_max, a_step, b_min, b_max, b_step, c_min, c_max, c_step);

    // gradient_descent_sine_wave(Io, n, &best_a, &best_b, &best_c, rate, iterations);

    double T = 2 * M_PI / best_b;

    printf("Hlavni poloosa: %f, Perioda: %f, Fazovy posun: %f\n", best_a, T, best_c);
    return 0;
}
