#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

// Параметры модели
const double a = 0.3;       // Коэффициент рождаемости жертв
const double b = 0.06;      // Коэффициент успешной охоты
const double gamma = 0.7;   // Коэффициент естественной убыли хищников
const double sigma = 0.035; // Коэффициент воспроизводства хищников

// Начальные условия
const double x0 = 5; // Начальное количество жертв
const double initial_y = 2; // Начальное количество хищников

// Функции для вычисления производных
double dx(double x, double y) {
    return (a - b * y) * x;
}

double dy(double x, double y) {
    return (-gamma + sigma * x) * y;
}

// Метод Рунге-Кутты 4-го порядка
void runge_kutta_method(double h, int steps, std::vector<double>& x, std::vector<double>& y) {
    for (int i = 0; i < steps; ++i) {
        double k1_x = h * dx(x[i], y[i]);
        double k1_y = h * dy(x[i], y[i]);

        double k2_x = h * dx(x[i] + k1_x / 2, y[i] + k1_y / 2);
        double k2_y = h * dy(x[i] + k1_x / 2, y[i] + k1_y / 2);

        double k3_x = h * dx(x[i] + k2_x / 2, y[i] + k2_y / 2);
        double k3_y = h * dy(x[i] + k2_x / 2, y[i] + k2_y / 2);

        double k4_x = h * dx(x[i] + k3_x, y[i] + k3_y);
        double k4_y = h * dy(x[i] + k3_x, y[i] + k3_y);

        double x_new = x[i] + (k1_x + 2 * k2_x + 2 * k3_x + k4_x) / 6;
        double y_new = y[i] + (k1_y + 2 * k2_y + 2 * k3_y + k4_y) / 6;

        x.push_back(x_new);
        y.push_back(y_new);
    }
}

int main() {
    // Параметры моделирования
    double h = 1.0; // Шаг дискретизации (1 день)
    int days = 365;  // Количество дней
    int steps = days; // Количество шагов моделирования

    // Векторы для хранения численности жертв и хищников
    std::vector<double> x = { x0 };
    std::vector<double> y = { initial_y };

    // Моделирование
    runge_kutta_method(h, steps, x, y);

    // Сохранение результатов для динамики жертв
    std::ofstream file_prey("prey_population.csv");
    file_prey << "Day,Prey\n";
    for (int i = 0; i <= steps; ++i) {
        file_prey << i << "," << x[i] << "\n";
    }
    file_prey.close();

    // Сохранение результатов для динамики хищников
    std::ofstream file_predator("predator_population.csv");
    file_predator << "Day,Predator\n";
    for (int i = 0; i <= steps; ++i) {
        file_predator << i << "," << y[i] << "\n";
    }
    file_predator.close();

    // Сохранение результатов для фазового портрета
    std::ofstream file_phase("phase_portrait.csv");
    file_phase << "Prey,Predator\n";
    for (int i = 0; i <= steps; ++i) {
        file_phase << x[i] << "," << y[i] << "\n";
    }
    file_phase.close();

    std::cout << "Моделирование завершено. Результаты сохранены в файлы:\n";
    std::cout << "1. prey_population.csv - динамика жертв\n";
    std::cout << "2. predator_population.csv - динамика хищников\n";
    std::cout << "3. phase_portrait.csv - фазовый портрет взаимодействия\n";

    return 0;
}
