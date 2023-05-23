/*
 * Tests.hpp
 *
 *  Created on: 26 срхїз 2019
 *      Author: Alex
 */

#ifndef TESTS_HPP_
#define TESTS_HPP_

std::vector<std::string> ChooseKLongest(const std::vector<std::string>& a, const int K);

void TestMostLikelySCS(const int reps, const int originalLen, const int copiesNum, const double delProb,
		std::mt19937& generator);
//void TestMostLikelyKCS(const int reps, const int originalLen, const double delProb, std::mt19937& generator);
void TestMostLikelyFilteredSCS(const int reps, const int originalLen, const int copiesNum, const int k,
		const double delProb, std::mt19937& generator);
void TestMinEditDistanceFilteredSCS(const int reps, const int originalLen, const int copiesNum, const int k,
		const double delProb, std::mt19937& generator);
void TestAlgorithm2(const int reps, const int originalLen, const int copiesNum, const double delProb,
		std::mt19937& generator);
void TestAlgorithm3(const int reps, const int originalLen, const int copiesNum, const double delProb,
		std::mt19937& generator);
void TestEitan(const int reps, const int originalLen, const double delProb, std::mt19937& generator);
void TestEitanOutput(const int reps, const int originalLen, const int copiesNum, const double delProb,
		std::mt19937& generator,const int maxOutputNum);
void TestEitanFull(const int reps, const int originalLen, const int copiesNum, const double delProb,
		std::mt19937& generator);
void TestSCSEstimate(const int reps, const int originalLen, const double delProb, std::mt19937& generator);
void TestBigClusterAlgorithm(const int reps, const int originalLen, const int copiesNum, const double delProb, const int k,
		std::mt19937& generator);
void TestAlgorithmNSteps(const int reps, const int maxStageNum, const int originalLen, const int copiesNum,
		const double delProb, std::mt19937& generator);
void TestAlgorithmNStepsHard(const int hardNum, const int maxStageNum, const int originalLen, const int copiesNum,
		const double delProb, std::mt19937& generator);
void TestSCSLen(const int reps, const int originalLen, const int copiesNum, const double delProb,
		std::mt19937& generator);
void TestCountPairwiseCorrectSize(const int reps, const int originalLen, const int copiesNum, const double delProb,
		std::mt19937& generator);
void TestFrequencies(const int reps, const int originalLen, const int copiesNum, const double delProb,
		std::mt19937& generator);
void TestSCSNum(const int reps, const int originalLen, const int K, const double delProb, std::mt19937& generator);
void TestKLongestSCSNum(const int reps, const int originalLen, const int copiesNum, const int K, const double delProb,
		std::mt19937& generator);
void TestMinSCSNum(const int reps, const int originalLen, const int copiesNum, const double delProb,
		std::mt19937& generator);
void TestChoseFromCorrectSizeSCSs(const int reps, const int originalLen, const int copiesNum, const double delProb,
		std::mt19937& generator);
void TestCorrectSizeFraction(const int reps, const int originalLen, const int copiesNum, const int k,
		const double delProb, std::mt19937& generator);
// KCS correctness tests

//void TestKCSNewAndOld(const int reps, const int minStrLen, const int maxStrLen, std::mt19937& generator);
//void TestKCSOriginalFoundFraction(const int reps, const int originalLen, const double delProb, std::mt19937& generator);
//void TestKCSIsSuperseq(const int reps, std::mt19937& generator, const int minStrLen, const int maxStrLen);
//void TestKCSLCS(const int reps, std::mt19937& generator, const int minStrLen, const int maxStrLen);
//void TestAllKCS(const int reps, std::mt19937& generator, const int minStrLen, const int maxStrLen);

// SCSN correctness tests
void Test1(int reps, std::mt19937& generator, const int minStrLen, const int maxStrLen);
void Test2(int reps, std::mt19937& generator, const int minStrLen, const int maxStrLen);
void TestSCSNIsSuperseq(int reps, std::mt19937& generator, const int minStrNum, const int maxStrNum,
		const int minStrLen, const int maxStrLen);
void TestSCSNInputOrder(int reps, std::mt19937& generator, const int minStrNum, const int maxStrNum,
		const int minStrLen, const int maxStrLen);
void TestSCSNHasAll(int reps, std::mt19937& generator, const int minStrNum, const int maxStrNum, const int minStrLen,
		const int maxStrLen, const int maxSCSLenToProcess);
void TestSCSNIsShortest(int reps, std::mt19937& generator, const int minStrNum, const int maxStrNum,
		const int minStrLen, const int maxStrLen, const int maxSCSLenToProcess);

void TestMostLikelySCSMaxEditDistance(const int reps, const int originalLen, const int copiesNum, const double delProb,
		std::mt19937& generator);
void TestMostLikelySCSKLongest(const int reps, const int originalLen, const int copiesNum, const int K,
		const double delProb, std::mt19937& generator);

// Timing tests

void TimeSCS2Len(const int reps, const int originalLen, const double delProb, std::mt19937& generator);
void TimeSCS2(const int reps, const int originalLen, const double delProb, std::mt19937& generator);
void TimeSCS2Multi(const int reps, const int originalLen, const double delProb, std::mt19937& generator);

#endif /* TESTS_HPP_ */
