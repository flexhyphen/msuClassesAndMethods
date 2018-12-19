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
    int numberOfPoint; // nomer tochki
    double x; // koordinata h
    double y; // koordinat u
    double z; // koordinata z
    double delta; // shum oblaka

    Point (double x_x = 0, double y_y = 0, double z_z = 0) { // konstruktor
        x = x_x;
        y = y_y;
        z = z_z;
    }

    Point (const Point &p) { // konstruktor
        x = p.x;
        y = p.y;
        z = p.z;
    }

    void operator = (const Point &point) { // operator prisvaivaniya
        x = point.x;
        y = point.y;
        z = point.z;
    }

    void operator += (const Point &point) { // operator +=
        x += point.x;
        y += point.y;
        z += point.z;
    }

    void operator -= (const Point &point) { // operator +=
        x -= point.x;
        y -= point.y;
        z -= point.z;
    }

    void operator *= (const Point &point) { //operator *=
        x *= point.x;
        y *= point.y;
        z *= point.z;
    }

    void operator /= (const Point &point) { //operator *=
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
    int n; // chislo tochek v oblake
    int numOfCloud; // personal'nyj nomer oblaka
    Point ctr; // centr
    //double xDsp, yDsp;
    Point DSP; // dispersiya
    Point *arrayOfPoints; // massiv tochek

    Cloud (int nn = 100, Point ctr_ctr = (static_cast<void>(0), static_cast<void>(0), 0), Point DSP_DSP = (static_cast<void>(0), static_cast<void>(0), 0) ) { // konstruktor
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


    void cloudStretching (int n, double e1, double e2) { // szhatie oblaka
        //cout << endl << "WORKED" << endl;
        for (int i = 0; i < n; i++) {
            arrayOfPoints[i] += Point(-ctr.x, -ctr.y);
            arrayOfPoints[i] *= Point(e1, e2);
            arrayOfPoints[i] += Point(ctr.x, ctr.y);
            //cout << endl << "WORKED" << i<< endl;
        }
    }

    void cloudMove (int n, Point point) { // sdvig oblaka na vektor
        for (int i = 0; i < n ; i++) {
            arrayOfPoints[i] += point;
        }
    }

    void cloudRotate (int n, double angle) { // povorot oblaka
        //cout << endl << angle;
        //cout << endl << M_PI;
        angle = angle * M_PI / 180;
        //cout << endl <<angle;
        double temp; // vspomogatel'naya peremennaya dlya zapominaniya predidushchego znacheniya
        for (int i = 0; i < n; i++) {
            temp = arrayOfPoints[i].x;
            arrayOfPoints[i].x = cos(angle) * temp - sin(angle) * arrayOfPoints[i].y;
            arrayOfPoints[i].y = sin(angle) * temp + cos(angle) * arrayOfPoints[i].y;
        }
    }

    void printCloud() { // pechat' oblaka v konsol'
        int j;
        // cout << endl << "N = " << n << endl; // otladka
        for (j = 0;j < n; j++) {
            cout << arrayOfPoints[j].x << " " << arrayOfPoints[j].y << " " << arrayOfPoints[j].z << endl;
        }
        //cout << endl << "WORKED2"<< endl; //otladka
    }
    void assignNumb() { // ustanovka nomera oblaka
        static int k = 0;
        k++;
        numOfCloud = k;
    }

    void printCloudInFile(int n, const char *fileName) { //pechat' oblaka v fajl
        ofstream fout(fileName);
        for (int k = 0; k < n; k++) {
            fout << arrayOfPoints[k].x << " " << arrayOfPoints[k].y << " " << arrayOfPoints[k].z << endl;
        }
    }

    Point centerOfMass(int n) { // poisk centra mass oblaka
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

    void cloudRotateCenterMass (int n, double angle) { // povorot otnositel'no centra mass
        //cout << endl << angle;
        //cout << endl << M_PI;
        Point pointCM = centerOfMass(n);
        angle = angle * M_PI / 180;
        //cout << endl <<angle;
        double temp;//  vspomogatel'naya peremennaya dlya zapominaniya predydushchego znacheniya
        for (int i = 0; i < n; i++) {
            temp = arrayOfPoints[i].x;
            arrayOfPoints[i].x = pointCM.x + cos(angle) * (temp - pointCM.x) - sin(angle) * (arrayOfPoints[i].y - pointCM.y);
            arrayOfPoints[i].y = pointCM.y + sin(angle) * (temp - pointCM.x) + cos(angle) * (arrayOfPoints[i].y - pointCM.y);
        }
    }


    void operator = (const Cloud &cloud) { // operator prisvaivaniya dlya oblaka
        ctr = cloud.ctr; // prisvaivaem centr
        DSP = cloud.DSP; // prisvaivaem pispersiyu
        n = cloud.n; // prisvaivaem chislo tochek
        numOfCloud = cloud.numOfCloud; // prisvaivaem nomer oblaka
        arrayOfPoints = new Point [n];
        for(int i = 0; i < n; i++) {
            arrayOfPoints[i] = cloud.arrayOfPoints[i];// prisvaivaem massiv tochek oblaka
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
    double A, B, C, D;// koehffy ploskosti
    double eigenValue1, eigenValue2, eigenValue3; // sobstvennye chisla
    Point eigenVector1, eigenVector2, eigenVector3; // sobstvennye vektora
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
                cloud.arrayOfPoints[i].z = s;
            }
            double sumx = 0;
            for (int j = 0; j < 100; j++) {
                sumx += (double)rand()*(5 + 5) / RAND_MAX + (-5);
            }
            double delta = sumx / 100;
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

        for(int i = 0; i < 3; i++){
            for (int j = 0; j < 3; j++){
                for(int k = 0; k < amountOfPoints; k++){
                    matrix[i][j] += x[k][i] * x[k][j];
                }
            }
        }

        findEigenValues(); // nahodim sobstvennye znacheniya
        eigenVector1 = eigenVector(matrix, eigenValue1);
        eigenVector2 = eigenVector(matrix, eigenValue2);
        eigenVector3 = eigenVector(matrix, eigenValue3);
        ofstream fout1("v1.txt");
        fout1 << cloud.ctr.x << " " << cloud.ctr.y << " " << cloud.ctr.z << " ";
        fout1 << cloud.ctr.x + eigenVector1.x << " " << cloud.ctr.y + eigenVector1.y << " " << cloud.ctr.z + eigenVector1.z ;

        ofstream fout2("v2.txt");
        fout2 << cloud.ctr.x << " " << cloud.ctr.y << " " << cloud.ctr.z << " ";
        fout2 << cloud.ctr.x + eigenVector2.x << " " << cloud.ctr.y + eigenVector2.y << " " << cloud.ctr.z + eigenVector2.z;

        ofstream fout3("v3.txt");
        fout3 << cloud.ctr.x << " " << cloud.ctr.y << " " << cloud.ctr.z << " ";
        fout3 << cloud.ctr.x + eigenVector3.x << " " << cloud.ctr.y + eigenVector3.y << " " << cloud.ctr.z + eigenVector3.z;

    }

    void findEigenValues(){
        double trace = 0; // sled
        for (int i = 0; i < 3; i++) {
            trace += matrix[i][i]; // ishchem sled
        }
        double q = trace / 3; //vspomogatel'naya peremennaya q
        double p1 = pow(matrix[0][1], 2) + pow(matrix[0][2], 2) + pow(matrix[1][2], 2); //vspomogatel'naya peremennaya p1
        double p2 = pow((matrix[0][0] - q), 2) + pow((matrix[1][1] - q), 2) + pow((matrix[2][2] - q), 2) + 2 * p1;//vspomogatel'naya peremennaya p2
        double p = sqrt(p2 / 6);//vspomogatel'naya peremennaya p
        double matrixB[3][3]; // vspomogatel'naya matrica poluchennaya
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

class Field { // klass pole

public:
    int numOfClouds; // kolichestvo oblakov
    int numOfComponents; // kolichestvo komponent
    int numOfClasters; // kolichestvo klasterov
    int numOfClastKernel; // kolichestvo klasterov s yadrami
    Cloud *arrayOfClouds; // massiv oblakov
    Cloud *arrayOfComp; // massiv komponent
    Cloud *arrayOfClasters; // massiv klasterov
    Cloud *arrayOfClastKern; // yadra
    double **distanceMatrix; // matrica rasstoyanij
    int **binaryMatrix; // binarnaya matrica
    Point *cklast; // massiv centrov klasterov
    Point *ctrOfClastKern; // massiv centrov klasterov posle k means kernel trick
    Field (int n = 0, int m = 0, int d = 0, int l = 0) { // konstruktor polya
        arrayOfClouds = new Cloud[100]; // sozdaem massiv na 100 potencial'nyh oblakov
        arrayOfComp = new Cloud [100]; // sozdaem massiv na 100 potencial'nyh komponent
        arrayOfClasters = new Cloud [100];
        arrayOfClastKern = new Cloud[100];
        numOfClouds = n;
        numOfComponents = m;
        numOfClasters = d;
        numOfClastKernel = l;
    }

    void addCloud(Cloud &cloud) { // komanda polyu dobavit' oblako
        arrayOfClouds[cloud.numOfCloud - 1] = cloud;
    }

    void moveCloud(int cloudNum, Point p) { // komanda polyu sdvinut' odno iz ego oblakov na vektor
        arrayOfClouds[cloudNum - 1].cloudMove(arrayOfClouds[cloudNum - 1].n, p);
    }

    void rotateCloud(int cloudNum, double ang) { // komanda polyu povernut' odno iz ego oblakov
        arrayOfClouds[cloudNum - 1].cloudRotate(arrayOfClouds[cloudNum - 1].n, ang);
    }

    void stretchingCloud(int cloudNum, Point p) { // komanda polyu szhat' odno iz ego oblakov
        arrayOfClouds[cloudNum - 1].cloudStretching(arrayOfClouds[cloudNum - 1].n, p.x, p.y);
    }


    void printCloud(int cloudNum) { // komanda polyu napechatat' odno iz ego oblakov v konsol'
        arrayOfClouds[cloudNum - 1].printCloud();
        //       cout << endl << "WORKED1"<< endl;
    }

    void printCloudInFile(int cloudNum, const char *fileName) { // komanda polyu napechatat' odno iz ego oblakov
        arrayOfClouds[cloudNum - 1].printCloudInFile(arrayOfClouds[cloudNum - 1].n, fileName);
    }

    void printCompInFile(int cloudNum, const char *fileName) { // komanda polyu napechatat' odnu iz ego komponent
        arrayOfComp[cloudNum - 1].printCloudInFile(arrayOfComp[cloudNum - 1].n, fileName);
    }

    void printClastInFile(int cloudNum, const char *fileName) { // komanda polyu napechatat' odnu iz ego klaster
        arrayOfClasters[cloudNum - 1].printCloudInFile(arrayOfClasters[cloudNum - 1].n, fileName);
    }

    void printKernInFile(int cloudNum, const char *fileName) { // komanda polyu napechatat' odnu iz ego klaster
        arrayOfClastKern[cloudNum - 1].printCloudInFile(arrayOfClastKern[cloudNum - 1].n, fileName);
    }

    void rotateCloudCenterMass(int cloudNum, double ang) { // komanda polyu povernut' odno iz ego oblakov otositel'no centra mass
        arrayOfClouds[cloudNum - 1].cloudRotateCenterMass(arrayOfClouds[cloudNum - 1].n, ang);
    }

    double distance(Point one, Point two) { // funkciya poiska rasstoyaniya mezhdu dvumya tochkami polya
        double sumX, sumY;
        sumX = abs(one.x - two.x);
        sumY = abs(one.y - two.y);
        return sqrt(pow(sumX,2) + pow(sumY,2));
    }

    int amountOfPointsInField() { // funkciya podscheta tochek v pole
        int amountOfPoints = 0;
        for (int i = 0; i < numOfClouds; i++) {
            amountOfPoints += arrayOfClouds[i].n;
        }
        return amountOfPoints;
    }

    void allFieldFill() { // funkciya zapolneniya poslednego oblaka polya (zanosim tuda vse sushchestvuyushchie oblaka po poryadku)
        int n = amountOfPointsInField();
        Point newCenterOfMass;
        int count = 0; // chislo tochek v
        int j = 0; // schetchik po ot 0 do kolichestva tochek v pole
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

    void distanceMatrixFill() { // zapolnenie matricy rasstoyanij
        int amountOfPoints = amountOfPointsInField(); // chislo tochek polya
        allFieldFill(); // zapolnenie poslednego oblaka vsemi predidushchimi oblakami
        distanceMatrix = new double *[amountOfPoints]; // matrica rasstoyanij
        for (int i = 0; i < amountOfPoints; i++) {
            distanceMatrix[i] = new double [amountOfPoints];
        }
        for (int i = 0; i < amountOfPoints; i++){
            for (int k = i; k < amountOfPoints; k++) {
                distanceMatrix[k][i] = distanceMatrix[i][k] =  distance(arrayOfClouds[99].arrayOfPoints[i], arrayOfClouds[99].arrayOfPoints[k]);
                //cout << distanceMatrix[i][k] << "otladka" << endl;
            }
        }
    }

    void printDistanceMatrix() { // pechat' dvoichnoj matricy dlya otladki
        int amountOfPoints = amountOfPointsInField();
        cout << "Matrica rasstoyanij: " << endl;
        distanceMatrixFill();
        for (int i = 0; i < amountOfPoints; i++){
            for (int k = 0; k < amountOfPoints; k++) {
                cout << distanceMatrix[i][k] << " ";
            }
            cout << endl;
        }
    }

    void binaryMatrixFill(double threshold) { // zapolnenie binarnoj matricy
        int amountOfPoints = amountOfPointsInField(); // chislo tochek polya
        distanceMatrixFill(); // zapolnyaem matricu rasstoyanij
        binaryMatrix = new int *[amountOfPoints];
        for (int i = 0; i < amountOfPoints; i++) {
            binaryMatrix[i] = new int [amountOfPoints];
        }
        for (int i = 0; i < amountOfPoints; i++) {
            for (int k = i; k < amountOfPoints; k++) {
                if( k == i) {
                    binaryMatrix[i][k] = binaryMatrix[k][i] = 0; // 0 na diagonali
                }
                else if(distanceMatrix[i][k] < threshold){
                    binaryMatrix[i][k] = binaryMatrix[k][i] = 1; // 1 esli rasstoyanie men'she poroga
                }
                else {
                    binaryMatrix[i][k] = binaryMatrix[k][i] = 0; // 0 esli rasstoyanie bol'she poroga
                }
            }
        }
    }

    void waveAlgorithm(const Cloud &allField) { // volnovoj algoritm
        int amountOfPointsInAllField = amountOfPointsInField(); // kolichestvo tochek v pole
        int i; // vspomogatel'nyj schetchik
        int l = 0; // vspomogatel'nyj schetchik
        int counter; // special'nyj schetchik dlya soseda
        numOfComponents = 0; // obnulyaem chislo komponent
        int pIOC = 0; // schetchik tochek v odnoj komponente
        int passed = 0; // schetchik projdennyh algoritmom tochek
        int *waveMatrix = new int [amountOfPointsInAllField]; // sozdaem vektor dlya poiska komponent tochek
        for (i = 0; i < amountOfPointsInAllField; i++) {
            waveMatrix[i] = 0; // zapolnenie vektora nulyami tak kak 0 shagov projdeno ( nachal'noe polozhenie )
        }
        int step; // takt
        bool change = false; // flazhok
        while (passed < amountOfPointsInAllField) { // poka kolichestvo projdenyh tochek men'she chem chislo tochek v pole
            step = 1; // pervyj takt vsego pri t=1
            for (i = 0; i < amountOfPointsInAllField; i++) {
                if (waveMatrix[i] == 0) { // ishchem neprojdennuyu tochku v massive tochek
                    waveMatrix[i] = step; // zadaem shag
                    step++; // inkrementiruem takt
                    pIOC++; // inkrementiruem chislo tochek v komponente
                    break;
                }
            }
            for (counter = 0; counter < amountOfPointsInAllField; counter++) { // ishchem sosedej
                change = false;
                if (waveMatrix[counter] == step - 1) { // esli tochka sosed nachal'noj tochki
                    for (i = 0; i < amountOfPointsInAllField; i++) {
                        if (binaryMatrix[counter][i] == 1 && waveMatrix[i] == 0) { // proveryaem
                            waveMatrix[i] = step; // prisvaivaem sosedu nomer takta
                            pIOC++; // inkrementiruem chislo tochek v komponente
                            change = true; // sosedi najdeny znachit ok
                        }
                    }
                }
                if (change) step++; // esli najden hot' odin sosed
            }
            passed += pIOC; // uvelchivaem chislo projdennyh tochek
            Cloud cloudComp(pIOC, (static_cast<void>(0),0), (static_cast<void>(1),1)); // sozdaem oblako - komponentu
            for (i = 0; i < amountOfPointsInAllField; i++) {
                if (waveMatrix[i] > 0) { // esli tochka byla najdena v takte
                    cloudComp.arrayOfPoints[l] = allField.arrayOfPoints[i]; // zanosim tochki prinadlezhshchie najdenoj komponente v oblako - komponentu
                    l++;
                    waveMatrix[i] = -1; // pomechaem vse projdenye tochki vyborki kak neaktivnye dlya novyh komponent
                }
            }
            cloudComp.numOfCloud = numOfComponents - 1; // daem oblaku-komponente personal'nyj nomer
            arrayOfComp[numOfComponents] = cloudComp; // zapisyvaem oblako-komponentu v massiv komponent
            numOfComponents++; // uvelichivaem znachenie chisla komponent
            pIOC = 0; // obnulyaem schetchik tochek v komponente
            if (change) step++;
            l = 0; // obnulyaem vspomogatel'nyj schetchik
        }
    }


    void printBinaryMatrix(double threshold) { // pechat' binarnoj matricy dlya otladki
        int amountOfPoints = amountOfPointsInField();
        binaryMatrixFill(threshold); // zapolnenie binarnoj matricy
        cout << "Binarnaya matrica rasstoyanij: " << endl;
        for (int i = 0; i < amountOfPoints; i++){
            for (int k = 0; k < amountOfPoints; k++){
                cout << binaryMatrix[i][k] << " ";
            }
            cout << endl;
        }
    }

    void KMeans(int k, const Cloud &cloud) {
        allFieldFill(); // zapolnyaem poslednee oblako vsemi oblakami
        int i = 0, j, t, g; // schetchiki
        int amountOfPoints = amountOfPointsInField(); // chislo tochek v pole
        int temp; // peremennaya dlya metki i proverki na ee smenu
        int *pointsInClast = new int [k]; // massiv kolichestva chisla tochek v klastere
        int *marks = new int [amountOfPoints]; // massiv prinadlezhnosti tochek polya k klasteram
        double threshold, xCTR,yCTR; // vspomogatel'nye peremennye dlya rasstoyaniya/h-koordinaty centra mass/y-k.c.m.
        bool change = true; //flazhok
        Point *centres = new Point [k]; // centry klasterov (starye)
        Point *ctrOfMass = new Point [k]; // centry klasterov (novye)

        for (i = 0; i < amountOfPoints; i++) {
            marks[i] = (-1); // zapolnyaem massiv markerov na tochkah (dlya otmetki centrov)
        }

        for (i = 0; i < k; i++) {
            g = rand() % amountOfPoints;
            marks[g] = i; // vybiraem markery centrov iz massiva prinadlezhnosti tochek k klasteram cherez kazhdye (aop/k tochek)
            centres[i] = cloud.arrayOfPoints[g]; // prisvaivaem nachal'nye centry po po indeksam markerov
        }

        while (change == true) { // poka centry menyayutsya
            for (i = 0; i < k; i++) {
                ctrOfMass[i] = Point(0, 0); // zapolnyaem massiv centrov mass
                pointsInClast[i] = 0; // zapolnyaem massiv kolichestva tochek
            }

            for (i = 0; i < amountOfPoints; i++) { // prohodim po kazhdoj tochke i otnosim ee k klasteru
                change = false; // menyaem flazhok
                temp = marks[i]; // zapominaem nachal'nuyu metku tochki
                threshold = distance(cloud.arrayOfPoints[i], centres[0]); //rasstoyanie mezhdu centrom pervogo klastera i i-oj tochkoj polya
                marks[i] = 0;

                for (j = 1; j < k; j++) {
                    if (distance(cloud.arrayOfPoints[i], centres[j]) < threshold) { //esli rasstoyanie ot i-oj tochki do j-ogo centra men'she nachal'nogo
                        marks[i] = j; // i-aya tochka prinadlezhit j-omu centru
                        threshold = distance(cloud.arrayOfPoints[i], centres[j]); // menyaem porog rasstoyaniya
                    }
                }

                ctrOfMass[marks[i]] += cloud.arrayOfPoints[i]; // summa po tochkam klastera
                pointsInClast[marks[i]]++; //  inkrementiruem chislo tochek v mark[i] klastere

                if (marks[i] != temp) { // esli markery centrov menyayutsya to prohodim po algoritmu eshche
                    change = true; // menyaem flag
                }
            }

            for (j = 0; j < k; j++) { // obnovlyaem centry iskhodya iz raspredeleniya tochek po klasteram
                xCTR = ctrOfMass[j].x / pointsInClast[j];
                yCTR = ctrOfMass[j].y / pointsInClast[j];
                ctrOfMass[j] = Point(xCTR,yCTR);
                centres[j] = ctrOfMass[j];
            }
        }

        for (i = 0; i < k; i++) { // zapolnyaem massiv klasterov
            Cloud claster(pointsInClast[i], (static_cast<void>(centres[i].x), centres[i].y), (static_cast<void>(1),1)); // sozdaem oblako-klaster
            t = 0; // vspomogatel'naya peremennaya
            for (j = 0; j < amountOfPoints; j++) { // smotrim massiv prinadlezhnostej tochek klasteram i
                if (marks[j] == i) { // esli j-aya tochka prenadlezhit i-omu klasteru
                    claster.arrayOfPoints[t] = cloud.arrayOfPoints[j]; // prisvaivaem massivu
                    t++;
                }
            }
            arrayOfClasters[i] = claster;  // prisvaivaem massiv klasterov
        }

        cklast = new Point[k]; // inicializaciya massiva centrov klasterov

        for (i = 0; i < k; i++) {
            cklast[i] = centres[i]; // zapolnyaem massiv centrov
        }

        delete [] ctrOfMass; delete [] pointsInClast; delete [] centres; delete [] marks; // chistka pamyati
        numOfClasters = k; // prisvaivaem chislo K klasterov atributu polya
    }

    void KMeansKernel(int k, int m, const Cloud &cloud) {
        void(KMeans(k, cloud)); //  iznachal'no razbivaem nashe mnozhestvo na K klasterov
        int i, j, h, amountOfPoints = amountOfPointsInField(), mark1;
        Cloud *sgust = new Cloud[k];
        Point *centres = new Point[m * k]; // inicializiruem centry s yadrami
        int *pointInClast = new int [k];
        int *marks = new int [amountOfPoints]; // massiv prinadlezhnosti tochek polya k klasteram
        double distation, temp, xCTR, yCTR; // porog rasstoyaniya / peremennaya dlya svapa rasstoyaniya / centr mass H / centr mass U
        bool change = true;
        Point *ctrOfMass = new Point [k];//centry klasterov
        for (i = 0; i < k; i++) {
            sgust[i] = arrayOfClasters[i]; // zapolnyaem vspomogatel'nyj massiv klasterov
        }

        for (i = 0; i < amountOfPoints; i++) {
            marks[i] = (-1); // zapolnyaem massiv prindlezhnosti tochek klasteram
        }
        for (i = 0; i < k; i++) {
            ctrOfMass[i] = Point(0, 0); // massiv centrov
            pointInClast[i] = 0; // zapolnyaem massiv chisla tochek v klastere
        }
        while (change == true) { // poka centry izmenyayutsya
            for (i = 0; i < k; i++) {
                ctrOfMass[i] = Point(0, 0); // massiv centrov
                pointInClast[i] = 0; // zapolnyaem massiv chisla tochek v klastere
            }
            for (i = 0; i < k; i++) {
                void(KMeans(m, sgust[i])); // razbivaem kazhdyj klaster na m podklasterov
                for(j = 0; j < m; j++) {
                    centres[i * m + j] = Point(cklast[j].x, cklast[j].y); // zapmnili centry
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
                    if (temp < distation) {//esli rasstoyanie ot i-oj tochki men'she do drugogo centra men'she
                    marks[i] = j;
                    distation = temp;
                    }
                    temp = 0;
                }
                ctrOfMass[marks[i]] += cloud.arrayOfPoints[i];
                pointInClast[marks[i]]++;
                if (marks[i] != mark1) { // metka pomenyalas' prodolzhaem
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
                arrayOfClastKern[i] = ClastKern; // zapolnyaem massiv klasterov oblakami-lasterami
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
    int readCommand() { // funkciya interfejsa chtenie komandy s konsoli i ee obrabotka

        cout << "Vvedite komandu:" ;
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
            printf(" CREATE - sozdat' oblako \n STRETCH - rastyanut' oblako \n DELETE - udalit' oblako \n SAVE - sohranit' oblako \n ROTATE - povernut' oblako  \n MOVE - sdvinut' oblako \n SAVEINFILE - sohranit' v fajl \n ROTATECM - povorot otnositel'no centra mass \n PRINTDA - pechat' matricy rasstoyanij \n PRINTBA - pechat' matricy rasstoyanij binarnoj \n WAVE - volnovoj algoritm \n SAVECOMP - sohranit' vse komponenty v odin fajl \n KMEANS - K-Srednih \n SAVECLASTIN - sohranit' vse komponenty v odin fajl \n KERNEL - K-means s yadrami \n SAVEKERNIN - pechat' klastera s yadrom \n MYORNOT - svoj ili chuzhoj \n PROJ - проекция облака собственые вектора \n");
            return 30;
        }
        if (command == "EXIT") return -1;

        cout << "KOMANDA NE NAJDENA" << endl;
        cout << "Poprobujte eshche raz" << endl;
        return -2;
    }

    Point readCentralPoint() { // funkciya interfejsa funkciya chteniya koordinat centra
        double x,y;
        cout << endl << "Vvedite koordinaty centra : " << endl;
        cout << "x = " ;
        cin >> x;
        cout << endl << "y = " ;
        cin >> y;
        Point p(x,y);
        //cout << p.x << p.y;
        return p;
    }

    Point readDispersion() { // funkciya interfejsa chteniya h i u dispersij
        double x,y;
        cout << endl << "Vvedite znacheniya dispersii : " << endl;
        cout << "xDsp = " ;
        cin >> x;
        cout << endl << "yDsp = " ;
        cin >> y;
        Point p(x,y);
        return p;
    }

    int readNumberOfPoints() { // funkciya interfejsa chtenie kol-va tochek v sozdavaem oblake
        int n;
        cout << endl << "Vvedite kolichestvo tochek v oblake : " << endl;
        cout << "n = " ;
        cin >> n;
        //cout << n; //  otladka
        return n;
    }

    int readCloudNumber(int N) { // funkciya interfejsa chtenie nomera oblaka
        int numberOfCloud;
        cout << " Vvedite nomer oblaka: " << endl;
        cin >> numberOfCloud;
        if (numberOfCloud > N && numberOfCloud != 100) {
            cout << "Oblaka ne sushchestvuet " << numberOfCloud << " > "  << N << endl;
            return -1;
        }
        return numberOfCloud;
    }

    int readCompNumber(int N) { // funkciya interfejsa chtenie nomera komponenty
        int numberOfCloud;
        cout << " Vvedite nomer oblaka: " << endl;
        cin >> numberOfCloud;
        if (numberOfCloud > N) {
            cout << "Oblaka ne sushchestvuet " << numberOfCloud << " > "  << N << endl;
            return -1;
        }
        return numberOfCloud;
    }

    int readClastNumber(int N) { // funkciya interfejsa chtenie nomera komponenty
        int numberOfCloud;
        cout << " Vvedite nomer oblaka: " << endl;
        cin >> numberOfCloud;
        if (numberOfCloud > N) {
            cout << "Oblaka ne sushchestvuet " << numberOfCloud << " > "  << N << endl;
            return -1;
        }
        return numberOfCloud;
    }

    Point readMove() {  // funkciya interfejsa schityvanie koordinat sdviga
        double x,y;
        cout << "Vvedite koordinaty sdviga : ";
        cout << "x = " << endl;
        cin >> x;
        cout << "y = " << endl;
        cin >> y;
        Point p(x,y);
        return p;
    }

    Point readStretch() { // funkciya interfejsa koehfficienty szhatiya
        double x,y;
        cout << "Vvedite koehfficienty rastyazheniya : " << endl;
        cout << "x = " ;
        cin >> x;
        cout << endl << "y = " ;
        cin >> y;
        Point p(x,y);
        return p;
    }

    double readAngle() { // // funkciya interfejsa chtenie ugla povorota v gradusah
        double angle;
        cout << "Vvedite ugol povorota: ";
        cout << endl << "angle = " ;
        cin >> angle;
        return angle;
    }

    double readThreshold() { // funkciya schityvaniya poroga dlya zapolneniya binarnoj matricy
        double threshold;
        cout << "Vvedite porog: ";
        cout << endl << "threshold = " ;
        cin >> threshold;
        return threshold;
    }

    int readK() { // funkciya schityvaniya poroga dlya zapolneniya binarnoj matricy
        int k;
        cout << "Vvedite K: ";
        cout << endl << "K = " ;
        cin >> k;
        return k;
    }

    int readKernel() { // funkciya schityvaniya poroga dlya zapolneniya binarnoj matricy
        int m;
        cout << "Vvedite CHislo YAder: ";
        cout << endl << "M = " ;
        cin >> m;
        return m;
    }

    double readPlane(string nameOfCoef) { // funkciya schityvaniya poroga dlya zapolneniya binarnoj matricy
        double i;
        cout << "Vvedite koehfficient " << nameOfCoef << ": ";
        cout << endl << "coefficient " << nameOfCoef << " = " ;
        cin >> i;
        return i;
    }

    string readFileName() {  // funkciya interfejsa chteniya nazvaniya fajla
        string fin;
        cout << "Vvedite imya fajla: ";
        cin >> fin;
        return fin;
    }

    Point readPoint() { // funkciya interfejsa chteniya tochki dlya SvojIliCHuzhoj
        double x,y;
        cout << "Vvedite koordinaty tochki dlya poiska nuzhnogo klastera : " << endl;
        cout << "x = " ;
        cin >> x;
        cout << endl << "y = " ;
        cin >> y;
        Point p(x,y);
        return p;
    }

};



class Controller { //klass kontroller

public:
    Interface mainInterface; // sozdaem ob"ekt odinochku - osnovnoj interfejs
    void createCloud(Field &field) { // funkciya
        int numOfPoints;
        Point ctr = mainInterface.readCentralPoint();
        Point Dsp = mainInterface.readDispersion();
        numOfPoints = mainInterface.readNumberOfPoints();
        //cout << endl << "N = "<<  numOfPoints << endl; //otladka
        Cloud cloudOne(numOfPoints, ctr, Dsp);
        cloudOne.assignNumb();
        field.addCloud(cloudOne);
        field.numOfClouds++;
        //cout << endl << "N = "<<  cloudOne.n << endl; //otladka
    }

    void moveCloud(Field &field) { // komanda ot kontrollera interfejsu i polyu dlya sdviga oblaka
        Point mov = mainInterface.readMove(); // chitaem vektor sdviga
        int n = mainInterface.readCloudNumber(field.numOfClouds); // chitaem nomer oblaka dlya sdviga
        field.moveCloud(n, mov); // daem komandu polyu sdvinut'
    }

    void stretchCloud(Field &field) { // komanda ot kontrollera interfejsu i polyu dlya rastyazheniya oblaka
        Point stretch = mainInterface.readStretch(); // chitaem koehfficienty rastyazheniya
        int n = mainInterface.readCloudNumber(field.numOfClouds); // chitaem nomer oblaka dlya rastyazheniya
        field.stretchingCloud(n, stretch); // daem komandu polyu rastyanut'
    }

    void rotateCloud(Field &field) { // komanda ot kontrollera interfejsu i polyu dlya povorota oblaka
        double ang = mainInterface.readAngle(); // chitaem ugol povorota
        int n = mainInterface.readCloudNumber(field.numOfClouds); // chitaem nomer oblaka dlya povorota
        field.rotateCloud(n, ang); // daem komandu polyu povernut' oblako
    }

    void printCloud(Field &field) { // komanda ot kontrollera interfejsu i polyu dlya pechati oblaka v konsol'
        int n = mainInterface.readCloudNumber(field.numOfClouds); // chitaem nomer oblaka dlya pechati
        field.printCloud(n); // daem komandu polyu napechat'
    }

    void rotateCloudCenterMass(Field &field) { //komanda ot kontrollera interfejsu i polyu dlya povotora otnositel'no centra mass
        double ang = mainInterface.readAngle();
        int n = mainInterface.readCloudNumber(field.numOfClouds);
        field.rotateCloudCenterMass(n, ang);
    }

    void printDistanceMatrix(Field &field) { // komanda ot kontrollera interfejsu i polyu dlya pechati matricy rasstoyanij
        //cout << field.amountOfPointsInField() << endl; // otladka
        field.printDistanceMatrix();
    }

    void printBinaryMatrix(Field &field) { // komanda ot kontrollera interfejsu i polyu dlya pechati binarnoj matricy
        double threashold = mainInterface.readThreshold();
        //cout << field.amountOfPointsInField() << endl; // otladka
        field.printBinaryMatrix(threashold);
    }

    void waveAlgorithm(Field &field) { // komanda ot kontrollera polyu dlya volnovogo algoritma
        field.waveAlgorithm(field.arrayOfClouds[99]);
    }

    void KMeans(Field &field) { // komanda ot kontrollera interfejsu ip polyu dlya K - srednih
        int k = mainInterface.readK();
        field.KMeans(k, field.arrayOfClouds[99]);
    }

    void KMeansKernel(Field &field) { // komanda ot kontrollera interfejsu i polyu dlya K - meand kernel trick
        int k = mainInterface.readK();
        int m = mainInterface.readKernel();
        field.KMeansKernel(k, m, field.arrayOfClouds[99]);
    }

    void allFieldFill(Field &field) { // komanda ot kontrollera polyu dlya zaneseniya vsekh tochek polya v odnoj oblako
        field.allFieldFill();
    }

    void projectionCloud(Field &field, Plane &plane) { //komanda ot kontrollera interfejsu i ploskosti dlya proekcirovaniya
        plane.A = mainInterface.readPlane("A");
        plane.B = mainInterface.readPlane("B");
        plane.C = mainInterface.readPlane("C");
        plane.D = mainInterface.readPlane("d");
        int n = mainInterface.readCloudNumber(field.numOfClouds);
        plane.matrixOfCloud(field.arrayOfClouds[n - 1]);
    }

    void myOrNot(Field &field) { // komanda ot kontrollera interfejsu i polyu dlya otveta na vopros kakomu klasteru prinadlezhit tochka
        Point unidefinedPoint = mainInterface.readPoint();
        char fin[100];
        string filename = mainInterface.readFileName();
        strcpy(fin, filename.c_str());
        unidefinedPoint.printPointInFile(fin);
        field.myOrNot(unidefinedPoint);
    }



    void createCloudInFile(Field &field, Point ctr, Point Dsp, int numOfPoints) { //komanda ot kontrollera interfejsu  i polyu dlya sozdaniya oblaka
        //cout << endl << "N = "<<  numOfPoints << endl; //otladka
        Cloud cloudOne(numOfPoints, ctr, Dsp);
        cloudOne.assignNumb();
        field.addCloud(cloudOne);
        field.numOfClouds++;
        //cout << endl << "N = "<<  cloudOne.n << endl; //otladka
    }

    void moveCloudInFile(Field &field, Point move, int n) { // komanda ot kontrollera interfejsu  i polyu dlya sdviga oblaka
        field.moveCloud(n, move); // daem komandu polyu sdvinut'
    }

    void stretchCloudInFile(Field &field, Point stretch, int n) { // komanda ot kontrollera interfejsu  i polyu dlya rastyazheniya oblaka
        field.stretchingCloud(n, stretch); // daem komandu polyu rastyanut'
    }

    void rotateCloudInFile(Field &field, double angle, int n) { // komanda ot kontrollera interfejsu  i polyu dlya povorota oblaka
        field.rotateCloud(n, angle); // daem komandu polyu povernut' oblako
    }

    void printCloudInFile(Field &field, int n) { // // komanda ot kontrollera interfejsu i polyu dlya pechati oblaka v konsol'
        field.printCloud(n); // daem komandu polyu napechat'
    }

    void rotateCloudCenterMassInFile(Field &field, double angle, int n) { //komanda ot kontrollera interfejsu i polyu dlya povotora otnositel'no centra mass
        field.rotateCloudCenterMass(n, angle);
    }

    void printDistanceMatrixInFile(Field &field) {// komanda ot kontrollera interfejsu i polyu dlya pechati matricy rasstoyanij
        //cout << field.amountOfPointsInField() << endl; // otladka
        field.printDistanceMatrix();
    }

    void printBinaryMatrixInFile(Field &field, double threshold) { // komanda ot kontrollera interfejsu i polyu dlya pechati binarnoj matricy
        //cout << field.amountOfPointsInField() << endl; // otladka
        field.printBinaryMatrix(threshold);
    }

    void waveAlgorithmInFile(Field &field) { // komanda ot kontrollera interfejsu ip polyu volnovoj algoritm
        field.waveAlgorithm(field.arrayOfClouds[99]);
    }

    void KMeansInFile(Field &field, int k) { // komanda ot kontrollera interfejsu ip polyu dlya K - srednih
        field.KMeans(k, field.arrayOfClouds[99]);
    }

    void KMeansKernelInFile(Field &field, int k, int m) { // komanda ot kontrollera interfejsu i polyu dlya K - meand kernel trick
        field.KMeansKernel(k, m, field.arrayOfClouds[99]);
    }

    void allFieldFillInFile(Field &field) { // komanda ot kontrollera polyu dlya zaneseniya vsekh tochek polya v odno oblako
        field.allFieldFill();
    }

    void projectionCloudInFile(Field &field, Plane &plane, int n, int A, int B, int C, int D) { //komanda ot kontrollera interfejsu i ploskosti dlya proekcirovaniya
        plane.A = A;
        plane.B = B;
        plane.C = C;
        plane.D = D;
        plane.matrixOfCloud(field.arrayOfClouds[n - 1]);
    }

    void myOrNotInFile(Field &field) { // komanda ot kontrollera interfejsu i polyu dlya otveta na vopros kakomu klasteru prinadlezhit tochka
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
        Plane plane; // создаем объект класса плоскость
        double how; // как вводить команды
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
                        mainController.waveAlgorithmInFile(field);
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
                        int k;
                        fin >> k;
                        mainController.KMeansInFile(field, k);
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
                        int n;
                        double A, B, C, D;
                        fin >> n >> A >> B >> C >> D;
                        mainController.projectionCloudInFile(field, plane, n, A, B, C, D);
                    }
                    if (valueOfCommand == 17) { // к - средних с ядрами
                        int k, m;
                        fin >> k >> m;
                        mainController.KMeansKernelInFile(field, k, m);
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

