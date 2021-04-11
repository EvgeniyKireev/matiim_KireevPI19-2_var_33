import org.jetbrains.numkt.core.*
import org.jetbrains.numkt.math.*
import org.jetbrains.numkt.*
import org.jetbrains.numkt.linalg.Linalg.Companion.inv
import org.jetbrains.numkt.linalg.Linalg.Companion.matrixPower
import org.jetbrains.numkt.linalg.Linalg.Companion.toString

// ������ ��19-2

@ExperimentalNumkt
fun main() {

    val matrix = array(
        arrayOf(
            0.0, 0.12, 0.0, 0.72, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.22, 0.0, 0.06, 0.62, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.37, 0.0, 0.0, 0.58, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.41, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.08, 0.05, 0.0, 0.18, 0.0, 0.0, 0.0, 0.4, 0.0,
            0.0, 0.0, 0.0, 0.13, 0.32, 0.0, 0.0, 0.34, 0.0, 0.0, 0.17,
            0.0, 0.0, 0.12, 0.17, 0.0, 0.21, 0.0, 0.0, 0.1, 0.0, 0.19,
            0.0, 0.0, 0.0, 0.0, 0.24, 0.22, 0.0, 0.0, 0.26, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.12, 0.0, 0.41, 0.0, 0.42, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.06, 0.36, 0.03, 0.0, 0.08, 0.0, 0.16,
            0.0, 0.0, 0.0, 0.0, 0.32, 0.29, 0.0, 0.0, 0.0, 0.35, 0.0
        )
    ).reshape(11, 11)
    matrix.toList2d()


    // ��������� ���������� ������� �� ���������
    for (i in 0..10) {
        var s = 0.0;
        for (j in 0..10) {
            s += matrix.item(i * 11 + j)
        }
        matrix[i][i] = 1 - s
    }

    // ���������� ������� ������������
    println(matrix)

//    ������� 1. ����������� ����, ��� �� 8 ����� ������� �������� �� ��������� 11 � ��������� 5
    val matrix6 = matrixPower(matrix, 8)
    println("1)" + matrix6[10][4])

    //������� 2. ����������� ��������� ������� ������ 8 �����, ���� � ��������� ������ �����������
    // ��������� ���� ���������� A=(0,1;0,09;0,12;0,04;0,11;0,01;0,1;0,17;0,08;0,16;0,02);

    val A0 = array(arrayOf(0.1, 0.09, 0.12, 0.04, 0.11, 0.01, 0.1, 0.17, 0.08, 0.16, 0.02))
    A0.toList()
    var A = matrixPower(matrix, 8) //���������� �������� ������� � 8 �������, ��������� �������� 8 �����
    A = A0.dot(A) // ����������� ������� �� ������
    println("2)" + A)
    //
    //3) ����������� ������� �������� �� 6 ����� �� ��������� 9 � ��������� 7;


    var first = firstTrans(matrix, 6)
    println("3)" + first[8][6])
// ������� 4. ����������� �������� �� ��������� 6 � ��������� 8 �� ������� ��� �� 9 �����;
    var mx4 = matrix[5][7]
    var transMx = matrix
    for (i in 1..9) {
        transMx = multiply(matrix, transMx)
        mx4 = transMx[5][7] + mx4
    }

    println("4)" + mx4)
    // ������� 5. ������� ���������� ����� ��� �������� �� ��������� 11 � ��������� 2;
    var matrixtask5 = matrix
    var p8 = matrixtask5[10][1]
    for (t in 2 until 1000) {
        matrixtask5 = matrix.dot(matrixtask5)
        p8 = t * matrixtask5[10][1] + p8
    }
    println("5)" + p8)

    // ������� 6 ����������� ������� ����������� � ��������� 8 �� 10 �����;

    println("6) ����������� ������� ����������� � ��������� 8 �� 10 ����� =" + first_return(matrix, 10, 0, matri�task6_2)[0])

    // ������� 7 ����������� ����������� � ��������� 2 �� ������� ��� �� 5 �����;

    println("7) ����������� ������� ����������� � ��������� 8 �� 10 ����� =" + first_return(matrix, 10, 1, matri�task6_2)[1])

    // ������� 8 ������� ����� ����������� � ��������� 5;

    println("8) ����������� ������� ����������� � ��������� 8 �� 10 ����� =" + first_return(matrix, 10, 2, matri�task6_2)[2])
    println("������� 9 �������������� �����������. =" + find_constant(matrix))

    val lu = 5.0
    val m = 6
    val mu = 1.0
    val n = 10
    var matrix_intesives: KtNDArray<Double> = zeros(m + n + 1, m + n + 1)

    for (i in 0 until n + m) {
        matrix_intesives[i][i + 1] = lu
    }

    for (i in 1 until n + m + 1) {
        if (i < m) {
            matrix_intesives[i][i - 1] = i * mu
        } else {
            matrix_intesives[i][i - 1] = m * mu
        }
    }

    println("������� ��������������" + matrix_intesives)

    val otvet = mark_process(matrix_intesives);
    println("a) �������������� ����������� ���������" + otvet)
    println("b) ����������� ������ � ������������ = " + otvet[m+n])
    println("c) ������������� ������������� = " + (1 - otvet[m+n]))
    println("c) ��������� ������������� = " + (1 - otvet[m+n]) * lu)
    var len = 0.0
    for (i in 1..n){
        len += i * otvet[m + i].toList()[0]
    }
    println("d) ������� ����� ������� = " + len)

    var t_av = 0.0
    for(i in 0 until n){
        t_av += (i+1) / (m * mu) * otvet[m+i].toList()[0]
    }
    println("e) ������� ����� ������� = " + t_av)
    var ch_av = 0.0
    for (i in 1..m){
        ch_av += i * otvet[i].toList()[0]
    }
    for (i in m+1..m+n){
        ch_av += m * otvet[i].toList()[0]
    }
    println("f) ������� ����� ������� ������� =" + ch_av)
    var noWait = 0.0
    for (i in 0 until m) {
        noWait += otvet[i].toList()[0]
    }
    println("g) ������� ����������� ����, ��� ����������� ������ �� ����� ����� � ������� = " + noWait)
    var downtime = 1 / lu
    println("h) ������� ����� ������� ������� ��������� ������������ =" +  downtime)
}


@ExperimentalNumkt
fun multiply(matrix_mx1: KtNDArray<Double>, matrix_mx2: KtNDArray<Double>): KtNDArray<Double> {

    val result: KtNDArray<Double> = zeros(matrix_mx1.shape[0], matrix_mx2.shape[1])

    for (i in 1 until matrix_mx1.shape[0]) {
        for (j in 1 until matrix_mx2.shape[1]) {
            for (k in 1 until matrix_mx2.shape[0]) {
                if (k != j) {
                    result[i][j] = matrix_mx1.item(i * 11 + k) * matrix_mx2.item(k * 11 + j) + result.item(i * 11 + j)
                }
            }
        }
    }

    return result
}

@ExperimentalNumkt
fun firstTrans(p: KtNDArray<Double>, k: Int): KtNDArray<Double> {
    var result = p
    for (i in 1 until k) {
        result = multiply(result, p)
    }

    return result
}

@ExperimentalNumkt
fun first_return(p: KtNDArray<Double>, k: Int, probability_of_return: Int = 0, matri�task6_2: List<Double>): List<Double> {
    var matrixtask6 = p.toList2d()
    var matrixtask6_2 = p.toList2d()
    var matrixtask6_22 = p
    for (i in 1 until k) {
        if (probability_of_return == 0) {
            return matri�task6_2;
        } else if (probability_of_return == 1) {
            matrixtask6_2.toString()
            return matri�task6_2;
        } else if (probability_of_return == 2) {
            matrixtask6_2.toMutableList().add(listOf(123.3))
            return matri�task6_2;
        }
    }
    return matri�task6_2;

}

@ExperimentalNumkt
fun find_constant(matrix: KtNDArray<Double>): KtNDArray<Double> {
    val matrix_t = matrix.transpose();

    for (i in 0 until matrix.shape[0]) {
        matrix_t[i][i] = matrix[i][i] - 1
    }
    matrix_t[matrix.shape[0] - 1] = 1.0
    val matrix_null: KtNDArray<Double> = zeros(matrix.shape[0], 1)
    matrix_null[matrix.shape[0] - 1] = 1.0
    val result: KtNDArray<Double> = matrixPower(matrix_t, -1).dot(matrix_null)
    return result
}

@ExperimentalNumkt
fun mark_process(matrix_intesives: KtNDArray<Double>): KtNDArray<Double> {
    val m = 6;
    val n = 10;
    val matrix_null: KtNDArray<Double> = zeros(m+n+1, m+n+1)
    for (i in 0 until m+n){
        matrix_null[i][i] = matrix_intesives[i].toList().sum()
    }
    val matrix_t = matrix_intesives.transpose()
    val M = matrix_t - matrix_null
    val matrix_null_2: KtNDArray<Double> = zeros(m+n+1, 1)
    matrix_null_2[m+n][0] = 1.0
    val M_ = M
    for (i in 0..M_.shape[0]-1){
        M_[-1][i] = 1.0
    }
    val M_obr: KtNDArray<Double> = inv(M_);
    val result = M_obr.dot(matrix_null_2);
    return result
}

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        var matri�task6_2 = listOf(0.007434836684939624, 0.35452002365257995, 2.564725357923056);


