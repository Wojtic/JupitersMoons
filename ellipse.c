#include <stdio.h>
#include <math.h>
#include <float.h>

typedef struct
{
    double x;
    double y;
} Point;

double solve_kepler(double e, double M)
{
    double E = M;
    double tolerance = 1e-6;
    double E_new;

    while (1)
    {
        E_new = M + e * sin(E);
        if (fabs(E_new - E) < tolerance)
            break;
        E = E_new;
    }

    return E_new;
}

void rotate_point(double *x, double *y, double angle)
{
    double newX = *x * cos(angle) - *y * sin(angle);
    *y = *x * sin(angle) + *y * cos(angle);
    *x = newX;
}
/*
int main()
{
    double a, e, t, T, angle;
    double t_max, dt;

    printf("Enter semi-major axis (in multiples of Jupiters diameter): ");
    scanf("%lf", &a);
    printf("Enter eccentricity (0 <= e < 1): ");
    scanf("%lf", &e);
    printf("Enter period (in days): ");
    scanf("%lf", &T);
    printf("Enter angle (in radians): ");
    scanf("%lf", &angle);
    printf("Enter total time to simulate (t_max in days): ");
    scanf("%lf", &t_max);
    printf("Enter time step (dt in days): ");
    scanf("%lf", &dt);

    FILE *file = fopen("orbit_data.txt", "w");
    if (file == NULL)
    {
        printf("Error opening file!\n");
        return 1;
    }

    for (t = 0; t <= t_max; t += dt)
    {
        double M = (2 * M_PI / T) * t;

        double E = solve_kepler(e, M);

        double x = a * (cos(E) - e);
        double y = a * sqrt(1 - e * e) * sin(E);

        rotate_point(&x, &y, angle);

        fprintf(file, "%.6f %.6f\n", t, y);
    }

    fclose(file);

    printf("Orbit data saved to orbit_data.txt\n");

    return 0;
}*/

double model(double a, double e, double T, double angle, double t)
{
    double M = (2 * M_PI / T) * t;

    double E = solve_kepler(e, M);

    double x = a * (cos(E) - e);
    double y = a * sqrt(1 - e * e) * sin(E);

    rotate_point(&x, &y, angle);

    return x;
}

double mean_squared_error(Point points[], int n, double a, double e, double T, double dT, double angle)
{
    double error = 0.0;
    for (int i = 0; i < n; i++)
    {
        double diff = model(a, e, T, angle, points[i].x + dT) - points[i].y;
        error += diff * diff;
    }
    return error / n;
}

void iterate_ellipse(Point points[], int n, double *best_a, double *best_e,
                     double *best_T, double *best_dT, double *best_angle,
                     double a_min, double a_max, double a_step,
                     double e_min, double e_max, double e_step,
                     double T_min, double T_max, double T_step,
                     double dT_min, double dT_max, double dT_step,
                     double angle_min, double angle_max, double angle_step)
{
    double min_error = DBL_MAX;

    for (double a = a_min; a <= a_max; a += a_step)
    {
        printf("a: %f\n", a);
        for (double e = e_min; e <= e_max; e += e_step)
        {
            for (double T = T_min; T <= T_max; T += T_step)
            {
                for (double angle = angle_min; angle <= angle_max; angle += angle_step)
                {
                    for (double dT = 0; dT <= T; dT += dT_step)
                    {
                        double error = mean_squared_error(points, n, a, e, T, dT, angle);

                        if (error < min_error)
                        {
                            min_error = error;
                            *best_a = a;
                            *best_e = e;
                            *best_T = T;
                            *best_dT = dT;
                            *best_angle = angle;
                        }
                    }
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
    Point Callisto[] = {
        {0, 13.30},
        {0.9549, 13.50},
        {1.9514, 11.80},
        {2.9410, 8.49},
        {4.7535, 0.00},
        {8.9827, -13.40}};

    int n = sizeof(points) / sizeof(points[0]);

    double a_min = 13.2, a_max = 13.8, a_steps = 100, a_step = (a_max - a_min) / a_steps;
    double e_min = 0, e_max = 0.015, e_steps = 50, e_step = (e_max - e_min) / e_steps;
    double T_min = 16.5, T_max = 17, T_steps = 50, T_step = (T_max - T_min) / T_steps;
    double dT_min = 0, dT_max = 1, dT_steps = 10, dT_step = (dT_max - dT_min) / dT_steps;
    double angle_min = 0, angle_max = 2 * M_PI, angle_steps = 10, angle_step = (angle_max - angle_min) / angle_steps;

    double best_a = 0.0, best_e = 0.0, best_T = 0.0, best_angle = 0.0, best_dT = 0.0;

    iterate_ellipse(points, n, &best_a, &best_e, &best_T, &best_dT, &best_angle,
                    a_min, a_max, a_step, e_min, e_max, e_step, T_min, T_max, T_step,
                    dT_min, dT_max, dT_step,
                    angle_min, angle_max, angle_step);

    // double T = 2 * M_PI / best_b;

    printf("a: %f, T: %f, e: %f, dt: %f, uhel: %f\n", best_a, best_T, best_e, best_dT, best_angle);
    return 0;
}
