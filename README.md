Задача:
Написать сишный extension для питона, вычисляющий пересечение двух прямых.
На входе:
1) первая прямая, заданная двумя точками [a, b]
2) вторая прямая, заданная двумя точками [c, d]
3) флаг А
4) флаг Б
 
На выходе:
1) точка пересечения двух прямых (находится на середине кратчайшего отрезка, соединяющего эти две прямые)
2) Расстояние между прямыми
3) если флаг А true, то выдавать скрещиваются ли прямые, если рассматривать их в качестве лучей (луч начинается в первой точке (a или с) и идет в сторону второй (b или d)). 
Лучи скрещиваются - значит точка пересечения двух прямых, содержащих эти лучи (середина кратчайшего отрезка, соединяющего прямые из 1-го пункта) находится в передних полупространствах каждого луча.
4) если флаг Б true, то вычислять расстояние от результирующей точки до прямой, соединяющей точки a и c
5) если прямые параллельны, выдавать в качестве точки пересечения None. Остальной вывод не важен.
 
Пример:
import numpy as np
a = np.array([-42.3, 14.4, 10.7], np.float)
b = np.array([11.8, 37.6, 4.9], np.float)
c = np.array([56.4, -27.1, 7.3], np.float)
d = np.array([40.2, 15.5, 8.4], np.float)
 
crossing([a, b], [c, d], True, True)
>>> [array([28.746, 44.774, 6.137]), 6.072, True, 55.61]
 
crossing([a, b], [c, d], False, False)
>>> [array([28.746, 44.774, 6.137]), 6.072, None, None]
