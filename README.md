# BBC_Tripipsiny
 
Назначение файлов:
 SourceFiles/CleanCoordx.xlsx - массив координат CpG островков в первой хромосоме (начало /t конец) в формате Excel
 SourceFiles/cpgIslandExt_chr1.tsv - tsv таблица с координатами CpG островков в первой хромосоме. Основа CleanCoords
 SourceFiles/generated_sequence.fa - меньший FASTA файл. Данные для первой задачи
 SourceFiles/hg38_chr1_and_chr2.fa - большой FASTA файл. Данные о двух хромосомах для второй задачи
 SourceFiles/myFASTA.fa - тестовый файл FASTA. Используется для дебага
 cgAllocation.py - тестовый скрипт. Выводит график, показывающий "концентрацию" CpG островков на небольшом участке нуклеотидной последовательности. Используется для тестового поиска CpG островков
 MakingMMatrix.py - скрипт, вычисляющий матрицу переходов Марковской цепи для заданной цепочки нуклеотидов. Первая часть решения второй задачи
 MatrixComparison.py - скрипт - решение первой задачи
 NucleotidCounter.py - тестовый скрипт. Выводит график с распределением нуклеотидов по всей нуклеотидной цепи
