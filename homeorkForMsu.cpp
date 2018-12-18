
//
//  main.cpp
//
//  Created by Alexander Vozvyshaev on 05.05.2018.
//  Copyright © 2018 Alexander Vozvyshaev. All rights reserved.
//

#include <cmath>
#include <cstdio>
#include <vector>
#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <cstring>
using namespace std;

class Point {
    
public:
    int numberOfPoint; // номер точки
    double x; // координата х
    double y; // координат у
    double z; // координата z
    double delta; // шум облака
    
    Point (double x_x = 0, double y_y = 0, double z_z = 0) { // конструктор
        x = x_x;
        y = y_y;
        z = z_z;
    }
    
    Point (const Point &p) { // конструктор
        x = p.x;
        y = p.y;
        z = p.z;
    }
    
    void operator = (const Point &point) { // оператор присваивания
        x = point.x;
        y = point.y;
        z = point.z;
    }
    
    void operator += (const Point &point) { // оператор +=
        x += point.x;
        y += point.y;
        z += point.z;
    }
    
    void operator -= (const Point &point) { // оператор +=
        x -= point.x;
        y -= point.y;
        z -= point.z;
    }
    
    void operator *= (const Point &point) { //оператор *=
        x *= point.x;
        y *= point.y;
        z *= point.z;
    }
    
    void operator /= (const Point &point) { //оператор *=
        x /= point.x;
        y /= point.y;
        z /= point.z;
    }
    void printPointInFile(const char *filename) {
        ofstream fout(filename);
        fout << x << " " << y << " " << endl;
    }
};

class Cloud {
    
public:
    int n; // число точек в облаке
    int numOfCloud; // персональный номер облака
    Point ctr; // центр
    //double xDsp, yDsp;
    Point DSP; // дисперсия
    Point *arrayOfPoints; // массив точек
    
    Cloud (int nn = 100, Point ctr_ctr = (static_cast<void>(0), static_cast<void>(0), 0), Point DSP_DSP = (static_cast<void>(0), static_cast<void>(0), 0) ) { // конструктор
        int k = 0;
        double x = 0, y = 0, z = 0;
        numOfCloud = k;
        n = nn;
        ctr = ctr_ctr;
        DSP = DSP_DSP;
        arrayOfPoints = new Point[nn];
        for(int j = 0; j < n; j++) {
            for (int i = 0; i < 1000; i++) {
                x += (double)rand() * (ctr.x + DSP.x - ctr.x + DSP.x) / RAND_MAX + ctr.x - DSP.x;
                y += (double)rand() * (ctr.y + DSP.y - ctr.y + DSP.y) / RAND_MAX + ctr.y - DSP.y;
            }
            x /= 1000;
            y /= 1000;
            arrayOfPoints[j] = Point(x, y, z);
            x = 0;
            y = 0;
        }
    }
    
    
    void cloudStretching (int n, double e1, double e2) { // сжатие облака
        //cout << endl << "WORKED" << endl;
        for (int i = 0; i < n; i++) {
            arrayOfPoints[i] += Point(-ctr.x, -ctr.y);
            arrayOfPoints[i] *= Point(e1, e2);
            arrayOfPoints[i] += Point(ctr.x, ctr.y);
            //cout << endl << "WORKED" << i<< endl;
        }
    }
    
    void cloudMove (int n, Point point) { // сдвиг облака на вектор
        for (int i = 0; i < n ; i++) {
            arrayOfPoints[i] += point;
        }
    }
    
    void cloudRotate (int n, double angle) { // поворот облака
        //cout << endl << angle;
        //cout << endl << M_PI;
        angle = angle * M_PI / 180;
        //cout << endl <<angle;
        double temp; // вспомогательная переменная для запоминания предидущего значения
        for (int i = 0; i < n; i++) {
            temp = arrayOfPoints[i].x;
            arrayOfPoints[i].x = cos(angle) * temp - sin(angle) * arrayOfPoints[i].y;
            arrayOfPoints[i].y = sin(angle) * temp + cos(angle) * arrayOfPoints[i].y;
        }
    }
    
    void printCloud() { // печать облака в консоль
        int j;
        // cout << endl << "N = " << n << endl; // отладка
        for (j = 0;j < n; j++) {
            cout << arrayOfPoints[j].x << " " << arrayOfPoints[j].y << " " << arrayOfPoints[j].z << endl;
        }
        //cout << endl << "WORKED2"<< endl; //отладка
    }
    void assignNumb() { // установка номера облака
        static int k = 0;
        k++;
        numOfCloud = k;
    }
    
    void printCloudInFile(int n, const char *fileName) { //печать облака в файл
        ofstream fout(fileName);
        for (int k = 0; k < n; k++) {
            fout << arrayOfPoints[k].x << " " << arrayOfPoints[k].y << " " << arrayOfPoints[k].z << endl;
        }
    }
    
    Point centerOfMass(int n) { // поиск центра масс облака
        double sumX = 0;
        double sumY = 0;
        for (int i = 0; i < n; i++) {
            sumX += arrayOfPoints[i].x;
            sumY += arrayOfPoints[i].y;
        }
        sumX /= n;
        sumY /= n;
        Point point(sumX, sumY);
        return point;
    }
    
    void cloudRotateCenterMass (int n, double angle) { // поворот относительно центра масс
        //cout << endl << angle;
        //cout << endl << M_PI;
        Point pointCM = centerOfMass(n);
        angle = angle * M_PI / 180;
        //cout << endl <<angle;
        double temp;//  вспомогательная переменная для запоминания предыдущего значения
        for (int i = 0; i < n; i++) {
            temp = arrayOfPoints[i].x;
            arrayOfPoints[i].x = pointCM.x + cos(angle) * (temp - pointCM.x) - sin(angle) * (arrayOfPoints[i].y - pointCM.y);
            arrayOfPoints[i].y = pointCM.y + sin(angle) * (temp - pointCM.x) + cos(angle) * (arrayOfPoints[i].y - pointCM.y);
        }
    }
    

    void operator = (const Cloud &cloud) { // оператор присваивания для облака
        ctr = cloud.ctr; // присваиваем центр
        DSP = cloud.DSP; // присваиваем писперсию
        n = cloud.n; // присваиваем число точек
        numOfCloud = cloud.numOfCloud; // присваиваем номер облака
        arrayOfPoints = new Point [n];
        for(int i = 0; i < n; i++) {
            arrayOfPoints[i] = cloud.arrayOfPoints[i];// присваиваем массив точек облака
        }
    }
    
    void printCloudInFile3DwithNoise(const char *FileName) {
        ofstream fout(FileName);
        for (int k = 0; k < n; k++) {
            fout<<arrayOfPoints[k].x<<"   "<< arrayOfPoints[k].y << "   "<< arrayOfPoints[k].z + arrayOfPoints[k].delta << "\n";
        }
    }
    
};

class Plane {
public:
    double A, B, C, D;// коэффы плоскости
    double eigenValue1, eigenValue2, eigenValue3; // собственные числа
    Point eigenVector1, eigenVector2, eigenVector3; // собственные вектора
    double matrix[3][3];
    
    Plane(double a = 1, double b = 1, double c = 1, double d = 1) {
        A = a;
        B = b;
        C = c;
        D = d;
    }
    
    void matrixOfCloud(Cloud &cloud) {
        int amountOfPoints = cloud.n;
        double x[amountOfPoints][3];
        
        for (int i = 0; i < amountOfPoints; i++) {
            Point p = cloud.arrayOfPoints[i];
            double s = -(A * p.x + B * p.y + D) / C;
            if (abs(C) > 0.000001) {
                s = cloud.arrayOfPoints[i].z;
            }
            double sumx = 0;
            for (int j = 0; j < 1000; j++) {
                sumx += (-1) + (rand() % 10000) * 0.0001;
            }
            double delta = sumx / 1000;
            cloud.arrayOfPoints[i].delta = delta;
        }
        cloud.printCloudInFile3DwithNoise("main.txt");
        for(int i = 0; i < 3; i++){
            for (int j = 0; j < 3; j++) {
                matrix[i][j] = 0;
            }
        }
        
        for(int i = 0; i < amountOfPoints; i++) {
            x[i][0] = cloud.arrayOfPoints[i].x;
            x[i][1] = cloud.arrayOfPoints[i].y;
            x[i][2] = cloud.arrayOfPoints[i].z;
        }
        
        for(int i = 0; i < amountOfPoints; i++){
            for (int j = 0; j < 3; j++){
                printf(" x[%i][%i] %f  ", i, j, x[i][j] );
            }
            cout << endl;
        }
        
        for(int i = 0; i < 3; i++){
            for (int j = 0; j < 3; j++){
                for(int k = 0; k < amountOfPoints; k++){
                    matrix[i][j] += x[k][i] * x[k][j];
                }
            }
        }
        
        findEigenValues(); // находим собственные значения
        eigenVector1 = eigenVector(matrix, eigenValue1);
        eigenVector2 = eigenVector(matrix, eigenValue2);
        eigenVector3 = eigenVector(matrix, eigenValue3);
        ofstream fout1("v1.txt"/*,ofstream::app*/);{
            fout1 << cloud.ctr.x <<"\t" << cloud.ctr.y << "\t" << cloud.ctr.z << endl;
            fout1 << cloud.ctr.x + eigenVector1.x << "\t" << cloud.ctr.y + eigenVector1.y << "\t" << cloud.ctr.z + eigenVector1.z << "\n";
        }
        ofstream fout2("v2.txt"/*,ofstream::app*/);{
            fout1 << cloud.ctr.x << "\t" << cloud.ctr.y << "\t" << cloud.ctr.z << endl;
            fout1 << cloud.ctr.x + eigenVector2.x << "\t" << cloud.ctr.y + eigenVector2.y << "\t" << cloud.ctr.z + eigenVector2.z << "\n";
        }
        ofstream fout3("v3.txt"/*,ofstream::app*/);{
            fout1 << cloud.ctr.x << "\t" << cloud.ctr.y << "\t" << cloud.ctr.z << endl;
            fout1 << cloud.ctr.x + eigenVector3.x << "\t" << cloud.ctr.y + eigenVector3.y << "\t" << cloud.ctr.z + eigenVector3.z << "\n";
        }
    }
    
    void findEigenValues(){
        double trace = 0; // след
        for (int i = 0; i < 3; i++) {
            trace += matrix[i][i]; // ищем след
        }
        double q = trace / 3; //вспомогательная переменная q
        double p1 = pow(matrix[0][1], 2) + pow(matrix[0][2], 2) + pow(matrix[1][2], 2); //вспомогательная переменная p1
        double p2 = pow((matrix[0][0] - q), 2) + pow((matrix[1][1] - q), 2) + pow((matrix[2][2] - q), 2) + 2 * p1;//вспомогательная переменная p2
        double p = sqrt(p2 / 6);//вспомогательная переменная p
        double matrixB[3][3]; // вспомогательная матрица полученная
        for(int i = 0; i < 3; i++){
            for (int j = 0; j < 3; j++){
                if(i == j) {
                    matrixB[i][i] = (matrix[i][i] - q) / p;
                }
                else {
                    matrixB[i][j] = matrix[i][j] / p;
                }
            }
        }
        double detB = findDeterminant3x3(matrixB);
        double r = detB / 2;
        double phi;
        if (r <= (-1)) {
            phi = M_PI / 3;
        }
        else if (r >= 1){
            phi = 0;
        }
        else {
            phi = acos(r) / 3;
        }
        
        eigenValue1 = q + 2 * p * cos(phi);
        eigenValue3 = q + 2 * p * cos(phi + (2 * M_PI / 3));
        eigenValue2 = 3 * q - eigenValue1 - eigenValue3;
    }
    
    double findDeterminant3x3(double a[3][3]) {
        double det = 0;
        det = a[0][0] * a[1][1] * a[2][2] + a[0][1] * a[1][2] * a[2][0] + a[0][2] * a[1][0] * a[2][1] - a[2][0] * a[1][1] * a[0][2] - a[1][0] * a[0][1] * a[2][2] - a[2][1] * a[1][2] * a[0][0];
        return det;
    }
    
    Point eigenVector(double matr[3][3], double eigenValue) {
        Point vector;
        vector.z = 1;
        double oper[3][3];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                oper[i][j] = matr[i][j];
                if (i == j) {
                    oper[i][j] -= eigenValue;
                }
            }
        }
        vector.y = (-oper[1][2] * oper[0][0] + oper[0][2] * oper[0][1]) / (oper[1][1] * oper[0][0] - oper[0][1] * oper[0][1]);
        vector.x = (-oper[0][2] - oper[0][1] * vector.y) / oper[0][0];
        return vector;
    }
        
    
};

class Field { // класс поле
    
public:
    int numOfClouds; // количество облаков
    int numOfComponents; // количество компонент
    int numOfClasters; // количество кластеров
    int numOfClastKernel; // количество кластеров с ядрами
    Cloud *arrayOfClouds; // массив облаков
    Cloud *arrayOfComp; // массив компонент
    Cloud *arrayOfClasters; // массив кластеров
    Cloud *arrayOfClastKern; // ядра
    double **distanceMatrix; // матрица расстояний
    int **binaryMatrix; // бинарная матрица
    Point *cklast; // массив центров кластеров
    Point *ctrOfClastKern; // массив центров кластеров после k means kernel trick
    Field (int n = 0, int m = 0, int d = 0, int l = 0) { // конструктор поля
        arrayOfClouds = new Cloud[100]; // создаем массив на 100 потенциальных облаков
        arrayOfComp = new Cloud [100]; // создаем массив на 100 потенциальных компонент
        arrayOfClasters = new Cloud [100];
        arrayOfClastKern = new Cloud[100];
        numOfClouds = n;
        numOfComponents = m;
        numOfClasters = d;
        numOfClastKernel = l;
    }
    
    void addCloud(Cloud &cloud) { // команда полю добавить облако
        arrayOfClouds[cloud.numOfCloud - 1] = cloud;
    }
    
    void moveCloud(int cloudNum, Point p) { // команда полю сдвинуть одно из его облаков на вектор
        arrayOfClouds[cloudNum - 1].cloudMove(arrayOfClouds[cloudNum - 1].n, p);
    }
    
    void rotateCloud(int cloudNum, double ang) { // команда полю повернуть одно из его облаков
        arrayOfClouds[cloudNum - 1].cloudRotate(arrayOfClouds[cloudNum - 1].n, ang);
    }
    
    void stretchingCloud(int cloudNum, Point p) { // команда полю сжать одно из его облаков
        arrayOfClouds[cloudNum - 1].cloudStretching(arrayOfClouds[cloudNum - 1].n, p.x, p.y);
    }
    
    
    void printCloud(int cloudNum) { // команда полю напечатать одно из его облаков в консоль
        arrayOfClouds[cloudNum - 1].printCloud();
        //       cout << endl << "WORKED1"<< endl;
    }
    
    void printCloudInFile(int cloudNum, const char *fileName) { // команда полю напечатать одно из его облаков
        arrayOfClouds[cloudNum - 1].printCloudInFile(arrayOfClouds[cloudNum - 1].n, fileName);
    }
    
    void printCompInFile(int cloudNum, const char *fileName) { // команда полю напечатать одну из его компонент
        arrayOfComp[cloudNum - 1].printCloudInFile(arrayOfComp[cloudNum - 1].n, fileName);
    }
    
    void printClastInFile(int cloudNum, const char *fileName) { // команда полю напечатать одну из его кластер
        arrayOfClasters[cloudNum - 1].printCloudInFile(arrayOfClasters[cloudNum - 1].n, fileName);
    }
    
    void printKernInFile(int cloudNum, const char *fileName) { // команда полю напечатать одну из его кластер
        arrayOfClastKern[cloudNum - 1].printCloudInFile(arrayOfClastKern[cloudNum - 1].n, fileName);
    }
    
    void rotateCloudCenterMass(int cloudNum, double ang) { // команда полю повернуть одно из его облаков отосительно центра масс
        arrayOfClouds[cloudNum - 1].cloudRotateCenterMass(arrayOfClouds[cloudNum - 1].n, ang);
    }
    
    double distance(Point one, Point two) { // функция поиска расстояния между двумя точками поля
        double sumX, sumY;
        sumX = abs(one.x - two.x);
        sumY = abs(one.y - two.y);
        return sqrt(pow(sumX,2) + pow(sumY,2));
    }
    
    int amountOfPointsInField() { // функция подсчета точек в поле
        int amountOfPoints = 0;
        for (int i = 0; i < numOfClouds; i++) {
            amountOfPoints += arrayOfClouds[i].n;
        }
        return amountOfPoints;
    }
    
    void allFieldFill() { // функция заполнения последнего облака поля (заносим туда все существующие облака по порядку)
        int n = amountOfPointsInField();
        Point newCenterOfMass;
        int count = 0; // число точек в
        int j = 0; // счетчик по от 0 до количества точек в поле
        Cloud cloudAllPoints(n, (static_cast<void>(0),0), (static_cast<void>(1),1));
        for(int i = 0; i < numOfClouds; i++){
            count = arrayOfClouds[i].n;
            for(int k = 0; k < count; k++){
                //cout << endl << &arrayOfClouds[i].arrayOfPoints[k] << endl;
                cloudAllPoints.arrayOfPoints[j] = arrayOfClouds[i].arrayOfPoints[k];
                j++;
            }
        }
        //cout << endl << "PRINT CLOUD ALL POINTS" << endl;
        //cloudAllPoints.printCloud();
        newCenterOfMass = cloudAllPoints.centerOfMass(cloudAllPoints.n);
        //cout << endl << newCenterOfMass.x << " " <<  newCenterOfMass.y << endl;
        cloudAllPoints.ctr = newCenterOfMass;
        arrayOfClouds[99] = cloudAllPoints;
    }
    
    void distanceMatrixFill() { // заполнение матрицы расстояний
        int amountOfPoints = amountOfPointsInField(); // число точек поля
        allFieldFill(); // заполнение последнего облака всеми предидущими облаками
        distanceMatrix = new double *[amountOfPoints]; // матрица расстояний
        for (int i = 0; i < amountOfPoints; i++) {
            distanceMatrix[i] = new double [amountOfPoints];
        }
        for (int i = 0; i < amountOfPoints; i++){
            for (int k = i; k < amountOfPoints; k++) {
                distanceMatrix[k][i] = distanceMatrix[i][k] =  distance(arrayOfClouds[99].arrayOfPoints[i], arrayOfClouds[99].arrayOfPoints[k]);
                //cout << distanceMatrix[i][k] << "отладка" << endl;
            }
        }
    }
    
    void printDistanceMatrix() { // печать двоичной матрицы для отладки
        int amountOfPoints = amountOfPointsInField();
        cout << "Матрица расстояний: " << endl;
        distanceMatrixFill();
        for (int i = 0; i < amountOfPoints; i++){
            for (int k = 0; k < amountOfPoints; k++) {
                cout << distanceMatrix[i][k] << " ";
            }
            cout << endl;
        }
    }
    
    void binaryMatrixFill(double threshold) { // заполнение бинарной матрицы
        int amountOfPoints = amountOfPointsInField(); // число точек поля
        distanceMatrixFill(); // заполняем матрицу расстояний
        binaryMatrix = new int *[amountOfPoints];
        for (int i = 0; i < amountOfPoints; i++) {
            binaryMatrix[i] = new int [amountOfPoints];
        }
        for (int i = 0; i < amountOfPoints; i++) {
            for (int k = i; k < amountOfPoints; k++) {
                if( k == i) {
                    binaryMatrix[i][k] = binaryMatrix[k][i] = 0; // 0 на диагонали
                }
                else if(distanceMatrix[i][k] < threshold){
                    binaryMatrix[i][k] = binaryMatrix[k][i] = 1; // 1 если расстояние меньше порога
                }
                else {
                    binaryMatrix[i][k] = binaryMatrix[k][i] = 0; // 0 если расстояние больше порога
                }
            }
        }
    }
    
    void waveAlgorithm(const Cloud &allField) { // волновой алгоритм
        int amountOfPointsInAllField = amountOfPointsInField(); // количество точек в поле
        int i; // вспомогательный счетчик
        int l = 0; // вспомогательный счетчик
        int counter; // специальный счетчик для соседа
        numOfComponents = 0; // обнуляем число компонент
        int pIOC = 0; // счетчик точек в одной компоненте
        int passed = 0; // счетчик пройденных алгоритмом точек
        int *waveMatrix = new int [amountOfPointsInAllField]; // создаем вектор для поиска компонент точек
        for (i = 0; i < amountOfPointsInAllField; i++) {
            waveMatrix[i] = 0; // заполнение вектора нулями так как 0 шагов пройдено ( начальное положение )
        }
        int step; // такт
        bool change = false; // флажок
        while (passed < amountOfPointsInAllField) { // пока количество пройденых точек меньше чем число точек в поле
            step = 1; // первый такт всего при т=1
            for (i = 0; i < amountOfPointsInAllField; i++) {
                if (waveMatrix[i] == 0) { // ищем непройденную точку в массиве точек
                    waveMatrix[i] = step; // задаем шаг
                    step++; // инкрементируем такт
                    pIOC++; // инкрементируем число точек в компоненте
                    break;
                }
            }
            for (counter = 0; counter < amountOfPointsInAllField; counter++) { // ищем соседей
                change = false;
                if (waveMatrix[counter] == step - 1) { // если точка сосед начальной точки
                    for (i = 0; i < amountOfPointsInAllField; i++) {
                        if (binaryMatrix[counter][i] == 1 && waveMatrix[i] == 0) { // проверяем
                            waveMatrix[i] = step; // присваиваем соседу номер такта
                            pIOC++; // инкрементируем число точек в компоненте
                            change = true; // соседи найдены значит ок
                        }
                    }
                }
                if (change) step++; // если найден хоть один сосед
            }
            passed += pIOC; // увелчиваем число пройденных точек
            Cloud cloudComp(pIOC, (static_cast<void>(0),0), (static_cast<void>(1),1)); // создаем облако - компоненту
            for (i = 0; i < amountOfPointsInAllField; i++) {
                if (waveMatrix[i] > 0) { // если точка была найдена в такте
                    cloudComp.arrayOfPoints[l] = allField.arrayOfPoints[i]; // заносим точки принадлежщие найденой компоненте в облако - компоненту
                    l++;
                    waveMatrix[i] = -1; // помечаем все пройденые точки выборки как неактивные для новых компонент
                }
            }
            cloudComp.numOfCloud = numOfComponents - 1; // даем облаку-компоненте персональный номер
            arrayOfComp[numOfComponents] = cloudComp; // записываем облако-компоненту в массив компонент
            numOfComponents++; // увеличиваем значение числа компонент
            pIOC = 0; // обнуляем счетчик точек в компоненте
            if (change) step++;
            l = 0; // обнуляем вспомогательный счетчик
        }
    }
    
    
    void printBinaryMatrix(double threshold) { // печать бинарной матрицы для отладки
        int amountOfPoints = amountOfPointsInField();
        binaryMatrixFill(threshold); // заполнение бинарной матрицы
        cout << "Бинарная матрица расстояний: " << endl;
        for (int i = 0; i < amountOfPoints; i++){
            for (int k = 0; k < amountOfPoints; k++){
                cout << binaryMatrix[i][k] << " ";
            }
            cout << endl;
        }
    }
    
    void KMeans(int k, const Cloud &cloud) {
        allFieldFill(); // заполняем последнее облако всеми облаками
        int i = 0, j, t, g; // счетчики
        int amountOfPoints = amountOfPointsInField(); // число точек в поле
        int temp; // переменная для метки и проверки на ее смену
        int *pointsInClast = new int [k]; // массив количества числа точек в кластере
        int *marks = new int [amountOfPoints]; // массив принадлежности точек поля к кластерам
        double threshold, xCTR,yCTR; // вспомогательные переменные для расстояния/х-координаты центра масс/y-к.ц.м.
        bool change = true; //флажок
        Point *centres = new Point [k]; // центры кластеров (старые)
        Point *ctrOfMass = new Point [k]; // центры кластеров (новые)
        
        for (i = 0; i < amountOfPoints; i++) {
            marks[i] = (-1); // заполняем массив маркеров на точках (для отметки центров)
        }
        
        for (i = 0; i < k; i++) {
            g = rand() % amountOfPoints;
            marks[g] = i; // выбираем маркеры центров из массива принадлежности точек к кластерам через каждые (aop/k точек)
            centres[i] = cloud.arrayOfPoints[g]; // присваиваем начальные центры по по индексам маркеров
        }

        while (change == true) { // пока центры меняются
            for (i = 0; i < k; i++) {
                ctrOfMass[i] = Point(0, 0); // заполняем массив центров масс
                pointsInClast[i] = 0; // заполняем массив количества точек
            }
            
            for (i = 0; i < amountOfPoints; i++) { // проходим по каждой точке и относим ее к кластеру
                change = false; // меняем флажок
                temp = marks[i]; // запоминаем начальную метку точки
                threshold = distance(cloud.arrayOfPoints[i], centres[0]); //расстояние между центром первого кластера и i-ой точкой поля
                marks[i] = 0;
                
                for (j = 1; j < k; j++) {
                    if (distance(cloud.arrayOfPoints[i], centres[j]) < threshold) { //если расстояние от i-ой точки до j-ого центра меньше начального
                        marks[i] = j; // i-ая точка принадлежит j-ому центру
                        threshold = distance(cloud.arrayOfPoints[i], centres[j]); // меняем порог расстояния
                    }
                }
                
                ctrOfMass[marks[i]] += cloud.arrayOfPoints[i]; // сумма по точкам кластера
                pointsInClast[marks[i]]++; //  инкрементируем число точек в mark[i] кластере
                
                if (marks[i] != temp) { // если маркеры центров меняются то проходим по алгоритму еще
                    change = true; // меняем флаг
                }
            }
            
            for (j = 0; j < k; j++) { // обновляем центры исходя из распределения точек по кластерам
                xCTR = ctrOfMass[j].x / pointsInClast[j];
                yCTR = ctrOfMass[j].y / pointsInClast[j];
                ctrOfMass[j] = Point(xCTR,yCTR);
                centres[j] = ctrOfMass[j];
            }
        }

        for (i = 0; i < k; i++) { // заполняем массив кластеров
            Cloud claster(pointsInClast[i], (static_cast<void>(centres[i].x), centres[i].y), (static_cast<void>(1),1)); // создаем облако-кластер
            t = 0; // вспомогательная переменная
            for (j = 0; j < amountOfPoints; j++) { // смотрим массив принадлежностей точек кластерам и
                if (marks[j] == i) { // если j-ая точка пренадлежит i-ому кластеру
                    claster.arrayOfPoints[t] = cloud.arrayOfPoints[j]; // присваиваем массиву
                    t++;
                }
            }
            arrayOfClasters[i] = claster;  // присваиваем массив кластеров
        }

        cklast = new Point[k]; // инициализация массива центров кластеров
        
        for (i = 0; i < k; i++) {
            cklast[i] = centres[i]; // заполняем массив центров
        }
        
        delete [] ctrOfMass; delete [] pointsInClast; delete [] centres; delete [] marks; // чистка памяти
        numOfClasters = k; // присваиваем число К кластеров атрибуту поля
    }
    
    void KMeansKernel(int k, int m, const Cloud &cloud) {
        void(KMeans(k, cloud)); //  изначально разбиваем наше множество на К кластеров
        int i, j, h, amountOfPoints = amountOfPointsInField(), mark1;
        Cloud *sgust = new Cloud[k];
        Point *centres = new Point[m * k]; // инициализируем центры с ядрами
        int *pointInClast = new int [k];
        int *marks = new int [amountOfPoints]; // массив принадлежности точек поля к кластерам
        double distation, temp, xCTR, yCTR; // порог расстояния / переменная для свапа расстояния / центр масс Х / центр масс У
        bool change = true;
        Point *ctrOfMass = new Point [k];//центры кластеров
        for (i = 0; i < k; i++) {
            sgust[i] = arrayOfClasters[i]; // заполняем вспомогательный массив кластеров
        }
        
        for (i = 0; i < amountOfPoints; i++) {
            marks[i] = (-1); // заполняем массив приндлежности точек кластерам
        }
        for (i = 0; i < k; i++) {
            ctrOfMass[i] = Point(0, 0); // массив центров
            pointInClast[i] = 0; // заполняем массив числа точек в кластере
        }
        while (change == true) { // пока центры изменяются
            for (i = 0; i < k; i++) {
                ctrOfMass[i] = Point(0, 0); // массив центров
                pointInClast[i] = 0; // заполняем массив числа точек в кластере
            }
            for (i = 0; i < k; i++) {
                void(KMeans(m, sgust[i])); // разбиваем каждый кластер на m подкластеров
                for(j = 0; j < m; j++) {
                    centres[i * m + j] = Point(cklast[j].x, cklast[j].y); // запмнили центры
                }
            }
            
            for (i = 0; i < amountOfPoints; i++) {
                change = false;
                mark1 = marks[i];
                distation = temp = 0;
                for (h = 0; h < m; h++) {
                    distation += distance(cloud.arrayOfPoints[i], (centres[h]));
                }
                marks[i] = 0;
                for (j = 1; j < k; j++) {
                    for (h = 0; h < m; h++) {
                        temp += distance(cloud.arrayOfPoints[i], (centres[j * m + h]));
                    }
                    if (temp < distation) {//если расстояние от i-ой точки меньше до другого центра меньше
                    marks[i] = j;
                    distation = temp;
                    }
                    temp = 0;
                }
                ctrOfMass[marks[i]] += cloud.arrayOfPoints[i];
                pointInClast[marks[i]]++;
                if (marks[i] != mark1) { // метка поменялась продолжаем
                    change = true;
                }
            }

            for (j = 0; j < k; j++) {
                xCTR = ctrOfMass[j].x / pointInClast[j];
                yCTR = ctrOfMass[j].y / pointInClast[j];
                ctrOfMass[j] = Point(xCTR, yCTR);
            }
           
            int t;
            for(i = 0; i < k; i++) {
                Cloud ClastKern(pointInClast[i], (static_cast<void>(ctrOfMass[i].x), ctrOfMass[i].y), (static_cast<void>(1), 1));
                t = 0;
                for(j = 0; j < amountOfPoints; j++) {
                    if (marks[j] == i) {
                        ClastKern.arrayOfPoints[t] = cloud.arrayOfPoints[j];
                        t++;
                    }
                }
                arrayOfClastKern[i] = ClastKern; // заполняем массив кластеров облаками-ластерами
            }
        }
       
        ctrOfClastKern = new Point[k*m];
        for (i = 0; i < k * m; i++) {
            ctrOfClastKern[i] = centres[i];
        }
        numOfClastKernel = k;
        delete [] ctrOfMass; delete [] pointInClast; delete [] centres; delete [] marks; delete [] sgust;
      
    }
    
    void myOrNot(Point &undefinded){
        
    }
    
    

    
};

class Interface {
    
public:
    string command;
    int readCommand() { // функция интерфейса чтение команды с консоли и ее обработка
        
        cout << "Введите команду:" ;
        cin >> command;
        if (command == "CREATE") return 1;
        if (command == "STRETCH") return 2;
        if (command == "SAVEINFILE") return 7;
        if (command == "DELETE") return 3;
        if (command == "SAVE") return 4;
        if (command == "ROTATE") return 5;
        if (command == "MOVE") return 6;
        if (command == "ROTATECM") return 8;
        if (command == "PROJ") return 16;
        if (command == "PRINTDA") return 9;
        if (command == "PRINTBA") return 10;
        if (command == "WAVE") return 11;
        if (command == "SAVECOMPIN") return 12;
        if (command == "KMEANS") return 13;
        if (command == "SAVECLASTIN") return 14;
        if (command == "FIELDFILL") return 15;
        if (command == "KERNEL") return 17;
        if (command == "SAVEKERNIN") return 18;
        if (command == "MYORNOT") return 19;

        if (command == "HELP") {
            printf(" CREATE - создать облако \n STRETCH - растянуть облако \n DELETE - удалить облако \n SAVE - сохранить облако \n ROTATE - повернуть облако  \n MOVE - сдвинуть облако \n SAVEINFILE - сохранить в файл \n ROTATECM - поворот относительно центра масс \n PRINTDA - печать матрицы расстояний \n PRINTBA - печать матрицы расстояний бинарной \n WAVE - волновой алгоритм \n SAVECOMP - сохранить все компоненты в один файл \n KMEANS - К-Средних \n SAVECLASTIN - сохранить все компоненты в один файл \n KERNEL - K-means с ядрами \n SAVEKERNIN - печать кластера с ядром \n MYORNOT - свой или чужой \n");
            return 7;
        }
        if (command == "EXIT") return -1;
        
        cout << "КОМАНДА НЕ НАЙДЕНА" << endl;
        cout << "Попробуйте еще раз" << endl;
        return -2;
    }
    
    Point readCentralPoint() { // функция интерфейса функция чтения координат центра
        double x,y;
        cout << endl << "Введите координаты центра : " << endl;
        cout << "x = " ;
        cin >> x;
        cout << endl << "y = " ;
        cin >> y;
        Point p(x,y);
        //cout << p.x << p.y;
        return p;
    }
    
    Point readDispersion() { // функция интерфейса чтения х и у дисперсий
        double x,y;
        cout << endl << "Введите значения дисперсии : " << endl;
        cout << "xDsp = " ;
        cin >> x;
        cout << endl << "yDsp = " ;
        cin >> y;
        Point p(x,y);
        return p;
    }
    
    int readNumberOfPoints() { // функция интерфейса чтение кол-ва точек в создаваем облаке
        int n;
        cout << endl << "Введите количество точек в облаке : " << endl;
        cout << "n = " ;
        cin >> n;
        //cout << n; //  отладка
        return n;
    }
    
    int readCloudNumber(int N) { // функция интерфейса чтение номера облака
        int numberOfCloud;
        cout << " Введите номер облака: " << endl;
        cin >> numberOfCloud;
        if (numberOfCloud > N && numberOfCloud != 100) {
            cout << "Облака не существует " << numberOfCloud << " > "  << N << endl;
            return -1;
        }
        return numberOfCloud;
    }
    
    int readCompNumber(int N) { // функция интерфейса чтение номера компоненты
        int numberOfCloud;
        cout << " Введите номер облака: " << endl;
        cin >> numberOfCloud;
        if (numberOfCloud > N) {
            cout << "Облака не существует " << numberOfCloud << " > "  << N << endl;
            return -1;
        }
        return numberOfCloud;
    }
    
    int readClastNumber(int N) { // функция интерфейса чтение номера компоненты
        int numberOfCloud;
        cout << " Введите номер облака: " << endl;
        cin >> numberOfCloud;
        if (numberOfCloud > N) {
            cout << "Облака не существует " << numberOfCloud << " > "  << N << endl;
            return -1;
        }
        return numberOfCloud;
    }
    
    Point readMove() {  // функция интерфейса считывание координат сдвига
        double x,y;
        cout << "Введите координаты сдвига : ";
        cout << "x = " << endl;
        cin >> x;
        cout << "y = " << endl;
        cin >> y;
        Point p(x,y);
        return p;
    }
    
    Point readStretch() { // функция интерфейса коэффициенты сжатия
        double x,y;
        cout << "Введите коэффициенты растяжения : " << endl;
        cout << "x = " ;
        cin >> x;
        cout << endl << "y = " ;
        cin >> y;
        Point p(x,y);
        return p;
    }
    
    double readAngle() { // // функция интерфейса чтение угла поворота в градусах
        double angle;
        cout << "Введите угол поворота: ";
        cout << endl << "angle = " ;
        cin >> angle;
        return angle;
    }
    
    double readThreshold() { // функция считывания порога для заполнения бинарной матрицы
        double threshold;
        cout << "Введите порог: ";
        cout << endl << "threshold = " ;
        cin >> threshold;
        return threshold;
    }
    
    int readK() { // функция считывания порога для заполнения бинарной матрицы
        int k;
        cout << "Введите K: ";
        cout << endl << "K = " ;
        cin >> k;
        return k;
    }
    
    int readKernel() { // функция считывания порога для заполнения бинарной матрицы
        int m;
        cout << "Введите Число Ядер: ";
        cout << endl << "M = " ;
        cin >> m;
        return m;
    }
    
    double readPlane(string nameOfCoef) { // функция считывания порога для заполнения бинарной матрицы
        double i;
        cout << "Введите коэффициент " << nameOfCoef << ": ";
        cout << endl << "coefficient " << nameOfCoef << " = " ;
        cin >> i;
        return i;
    }
    
    string readFileName() {  // функция интерфейса чтения названия файла
        string fin;
        cout << "Введите имя файла: ";
        cin >> fin;
        return fin;
    }
    
    Point readPoint() { // функция интерфейса чтения точки для СвойИлиЧужой
        double x,y;
        cout << "Введите координаты точки для поиска нужного кластера : " << endl;
        cout << "x = " ;
        cin >> x;
        cout << endl << "y = " ;
        cin >> y;
        Point p(x,y);
        return p;
    }
    
};



class Controller { //класс контроллер
    
public:
    Interface mainInterface; // создаем объект одиночку - основной интерфейс
    void createCloud(Field &field) { // функция
        int numOfPoints;
        Point ctr = mainInterface.readCentralPoint();
        Point Dsp = mainInterface.readDispersion();
        numOfPoints = mainInterface.readNumberOfPoints();
        //cout << endl << "N = "<<  numOfPoints << endl; //отладка
        Cloud cloudOne(numOfPoints, ctr, Dsp);
        cloudOne.assignNumb();
        field.addCloud(cloudOne);
        field.numOfClouds++;
        //cout << endl << "N = "<<  cloudOne.n << endl; //отладка
    }
    
    void moveCloud(Field &field) { // команда от контроллера интерфейсу  и полю для сдвига облака
        Point mov = mainInterface.readMove(); // читаем вектор сдвига
        int n = mainInterface.readCloudNumber(field.numOfClouds); // читаем номер облака для сдвига
        field.moveCloud(n, mov); // даем команду полю сдвинуть
    }
    
    void stretchCloud(Field &field) { // команда от контроллера интерфейсу  и полю для растяжения облака
        Point stretch = mainInterface.readStretch(); // читаем коэффициенты растяжения
        int n = mainInterface.readCloudNumber(field.numOfClouds); // читаем номер облака для растяжения
        field.stretchingCloud(n, stretch); // даем команду полю растянуть
    }
    
    void rotateCloud(Field &field) { // команда от контроллера интерфейсу  и полю для поворота облака
        double ang = mainInterface.readAngle(); // читаем угол поворота
        int n = mainInterface.readCloudNumber(field.numOfClouds); // читаем номер облака для поворота
        field.rotateCloud(n, ang); // даем команду полю повернуть облако
    }
    
    void printCloud(Field &field) { // команда от контроллера интерфейсу и полю для печати облака в консоль
        int n = mainInterface.readCloudNumber(field.numOfClouds); // читаем номер облака для печати
        field.printCloud(n); // даем команду полю напечать
    }
    
    void rotateCloudCenterMass(Field &field) { //команда от контроллера интерфейсу и полю для повотора относительно центра масс
        double ang = mainInterface.readAngle();
        int n = mainInterface.readCloudNumber(field.numOfClouds);
        field.rotateCloudCenterMass(n, ang);
    }
    
    void printDistanceMatrix(Field &field) { // команда от контроллера интерфейсу и полю для печати матрицы расстояний
        //cout << field.amountOfPointsInField() << endl; // отладка
        field.printDistanceMatrix();
    }
    
    void printBinaryMatrix(Field &field) { // команда от контроллера интерфейсу и полю для печати бинарной матрицы
        double threashold = mainInterface.readThreshold();
        //cout << field.amountOfPointsInField() << endl; // отладка
        field.printBinaryMatrix(threashold);
    }
    
    void waveAlgorithm(Field &field) { // команда от контроллера полю для волнового алгоритма
        field.waveAlgorithm(field.arrayOfClouds[99]);
    }
    
    void KMeans(Field &field) { // команда от контроллера интерфейсу ип полю для К - средних
        int k = mainInterface.readK();
        field.KMeans(k, field.arrayOfClouds[99]);
    }
    
    void KMeansKernel(Field &field) { // команда от контроллера интерфейсу и полю для К - meand kernel trick
        int k = mainInterface.readK();
        int m = mainInterface.readKernel();
        field.KMeansKernel(k, m, field.arrayOfClouds[99]);
    }
    
    void allFieldFill(Field &field) { // команда от контроллера полю для занесения всех точек поля в одной облако
        field.allFieldFill();
    }
    
    void projectionCloud(Field &field, Plane &plane) { //команда от контроллера интерфейсу и плоскости для проекцирования 
        plane.A = mainInterface.readPlane("A");
        plane.B = mainInterface.readPlane("B");
        plane.C = mainInterface.readPlane("C");
        plane.D = mainInterface.readPlane("d");
        int n = mainInterface.readCloudNumber(field.numOfClouds);
        plane.matrixOfCloud(field.arrayOfClouds[n - 1]);
    }
    
    
    
    void createCloudInFile(Field &field, Point ctr, Point Dsp, int numOfPoints) { // функция
        //cout << endl << "N = "<<  numOfPoints << endl; //отладка
        Cloud cloudOne(numOfPoints, ctr, Dsp);
        cloudOne.assignNumb();
        field.addCloud(cloudOne);
        field.numOfClouds++;
        //cout << endl << "N = "<<  cloudOne.n << endl; //отладка
    }
    
    void moveCloudInFile(Field &field, Point move, int n) { // команда от контроллера интерфейсу  и полю для сдвига облака
        field.moveCloud(n, move); // даем команду полю сдвинуть
    }
    
    void stretchCloudInFile(Field &field, Point stretch, int n) { // команда от контроллера интерфейсу  и полю для растяжения облака
        field.stretchingCloud(n, stretch); // даем команду полю растянуть
    }
    
    void rotateCloudInFile(Field &field, double angle, int n) { // команда от контроллера интерфейсу  и полю для поворота облака
        field.rotateCloud(n, angle); // даем команду полю повернуть облако
    }
    
    void printCloudInFile(Field &field, int n) { // // команда от контроллера интерфейсу  и полю для печати облака в консоль
        field.printCloud(n); // даем команду полю напечать
    }
    
    void rotateCloudCenterMassInFile(Field &field, double angle, int n) {
        field.rotateCloudCenterMass(n, angle);
    }
    
    void printDistanceMatrixInFile(Field &field) {
        //cout << field.amountOfPointsInField() << endl; // отладка
        field.printDistanceMatrix();
    }
    
    void printBinaryMatrixInFile(Field &field, double threshold) {
        //cout << field.amountOfPointsInField() << endl; // отладка
        field.printBinaryMatrix(threshold);
    }
    
    void waveAlgorithmInFile(Field &field) {
        field.waveAlgorithm(field.arrayOfClouds[99]);
    }
    
    void KMeansInFile(Field &field, int k) {
        field.KMeans(k, field.arrayOfClouds[99]);
    }
    
    void KMeansKernelInFile(Field &field, int k, int m) {
        field.KMeansKernel(k, m, field.arrayOfClouds[99]);
    }
    
    void allFieldFillInFile(Field &field) {
        field.allFieldFill();
    }
    
    void projectionCloudInFile(Field &field, Plane &plane, int n, int A, int B, int C, int D) {
        plane.A = A;
        plane.B = B;
        plane.C = C;
        plane.D = D;
        plane.matrixOfCloud(field.arrayOfClouds[n - 1]);
    }
    
    void myOrNot(Field &field) {
        Point unidefined = mainInterface.readPoint();
        char fin[100];
        string filename = mainInterface.readFileName();
        strcpy(fin, filename.c_str());
        unidefined.printPointInFile(fin);
        field.myOrNot(unidefined);
    }
    
    void start() { // стартовая команда для запуска приложения
        int valueOfCommand = 0;
        Field field(0); // создаем объект поле
        Plane plane;
        double how;
        Controller mainController; // создаем класс одиночку - основной контроллер
        cout << "Как будем вводить команды 0 - консоль , 1 - файл :" << endl;
        cin >> how;
        if(how == 0){
            while (valueOfCommand != -1){ // читаем команды с консоли пока не будет EXIT
                valueOfCommand = mainInterface.readCommand();
                if (valueOfCommand == 1) { // создать облако
                    mainController.createCloud(field);
                }
                if (valueOfCommand == 2) { // сжать облако
                    mainController.stretchCloud(field);
                }
                if (valueOfCommand == 4) { // распечатать облако
                    mainController.printCloud(field);
                }
                if (valueOfCommand == 5) { // повернуть облако
                    mainController.rotateCloud(field);
                }
                if (valueOfCommand == 6) { // сдвинуть облако
                    mainController.moveCloud(field);
                }
                if (valueOfCommand == 7) { // печать облака в файл
                    char fin[100];
                    int n = mainInterface.readCloudNumber(field.numOfClouds);
                    string filename = mainInterface.readFileName();
                    strcpy(fin, filename.c_str());
                    field.printCloudInFile(n, fin);
                }
                if (valueOfCommand == 8) { // повернуть относительно центра масс
                    mainController.rotateCloudCenterMass(field);
                }
                if (valueOfCommand == 9) { // распечатать  матрицу расстояний
                    mainController.printDistanceMatrix(field);
                }
                if (valueOfCommand == 10) { // печать бинарной матрицы
                    mainController.printBinaryMatrix(field);
                }
                if (valueOfCommand == 11) { // волновой алгоритм
                    mainController.waveAlgorithm(field);
                }
                if (valueOfCommand == 12) { // печать компоненты в файл
                    char fin[100];
                    int n = mainInterface.readCompNumber(field.numOfComponents);
                    string filename = mainInterface.readFileName();
                    strcpy(fin, filename.c_str());
                    field.printCompInFile(n, fin);
                }
                if (valueOfCommand == 13) { // к - средних
                    mainController.KMeans(field);
                }
                if (valueOfCommand == 14) { // печать кластера в файл
                    char fin[100];
                    int n = mainInterface.readClastNumber(field.numOfClasters);
                    string filename = mainInterface.readFileName();
                    strcpy(fin, filename.c_str());
                    field.printClastInFile(n, fin);
                }
                if (valueOfCommand == 15) { // занесение всех облаков в последнее облако и печать его в файл (комнда для теста)
                    mainController.allFieldFill(field);
                    char fin[100];
                    string filename = mainInterface.readFileName();
                    strcpy(fin, filename.c_str());
                    field.printCloudInFile(100, fin);
                }
                if (valueOfCommand == 16) { // проекцировать облако
                    mainController.projectionCloud(field, plane);
                }
                if (valueOfCommand == 17) { // к - средних с ядрами
                    mainController.KMeansKernel(field);
                }
                if (valueOfCommand == 18) { // печать кластера с ядрами в файл
                    char fin[100];
                    int n = mainInterface.readClastNumber(field.numOfClastKernel);
                    string filename = mainInterface.readFileName();
                    strcpy(fin, filename.c_str());
                    field.printKernInFile(n, fin);
                }
                if (valueOfCommand == 19) {
                    mainController.myOrNot(field);
                }

            }
        }
        else if (how == 1) {
            string name;
            cout << "Введите имя файла: ";
            cin >> name;
            ifstream fin;
            fin.open(name.c_str());
            while (valueOfCommand != -1) {
                while (valueOfCommand != -1){ // читаем команды с консоли пока не будет EXIT
                    fin >> valueOfCommand;
                    if (valueOfCommand == 1) { // создать облако
                        Point ctr;
                        Point Dsp;
                        int n;
                        fin >> ctr.x >> ctr.y >> Dsp.x >> Dsp.y >> n;
                        mainController.createCloudInFile(field, ctr, Dsp, n);
                    }
                    if (valueOfCommand == 2) { // сжать облако
                        Point stretch;
                        int n;
                        fin >> stretch.x >> stretch.y >> n;
                        mainController.stretchCloudInFile(field, stretch, n);
                    }
                    if (valueOfCommand == 4) { // распечатать облако
                        int n;
                        fin >> n;
                        mainController.printCloudInFile(field, n);
                    }
                    if (valueOfCommand == 5) { // повернуть облако
                        double angle;
                        int n;
                        fin >> angle >> n;
                        mainController.rotateCloudInFile(field, angle, n);
                    }
                    if (valueOfCommand == 6) { // сдвинуть облако
                        Point move;
                        int n;
                        fin >> move.x >> move.y >> n;
                        mainController.moveCloudInFile(field, move, n);
                    }
                    if (valueOfCommand == 7) { // печать облака в файл
                        char d[100];
                        int n;
                        string filename;
                        fin >> n >> filename;
                        strcpy(d, filename.c_str());
                        field.printCloudInFile(n, d);
                    }
                    if (valueOfCommand == 8) { // повернуть относительно центра масс
                        double angle;
                        int n;
                        fin >> angle >> n;
                        mainController.rotateCloudCenterMassInFile(field, angle, n);
                    }
                    if (valueOfCommand == 9) { // распечатать  матрицу расстояний
                        mainController.printDistanceMatrixInFile(field);
                    }
                    if (valueOfCommand == 10) { // печать бинарной матрицы
                        double threshold;
                        fin >> threshold;
                        mainController.printBinaryMatrixInFile(field, threshold);
                    }
                    if (valueOfCommand == 11) { // волновой алгоритм
                        mainController.waveAlgorithm(field);
                    }
                    if (valueOfCommand == 12) { // печать компоненты в файл
                        char d[100];
                        int n;
                        string filename;
                        fin >> n >> filename;
                        strcpy(d, filename.c_str());
                        field.printCompInFile(n, d);
                    }
                    if (valueOfCommand == 13) { // к - средних
                        mainController.KMeans(field);
                    }
                    if (valueOfCommand == 14) { // печать кластера в файл
                        char d[100];
                        int n;
                        string filename;
                        fin >> n >> filename;
                        strcpy(d, filename.c_str());
                        field.printClastInFile(n, d);
                    }
                    if (valueOfCommand == 15) { // занесение всех облаков в одно последнее облако и печать его (команда для теста)
                        mainController.allFieldFill(field);
                        char d[100];
                        string filename;
                        fin >> filename;
                        strcpy(d, filename.c_str());
                        field.printCloudInFile(100, d);
                    }
                    if (valueOfCommand == 16) { // проекцировать облако
                        mainController.projectionCloud(field, plane);
                    }
                    if (valueOfCommand == 17) { // к - средних с ядрами
                        mainController.KMeansKernel(field);
                    }
                    if (valueOfCommand == 18) { // печать кластера в файл
                        char d[100];
                        int n;
                        string filename;
                        fin >> n >> filename;
                        strcpy(d, filename.c_str());
                        field.printKernInFile(n, d);
                    }
                    if (valueOfCommand == 19) {
                        
                    }
                }
            }
            fin.close();
        }
    }
};

int main(void) {
    setlocale(LC_ALL,"Rus");
    Controller mainController;
    mainController.start();
    return 0;
}
