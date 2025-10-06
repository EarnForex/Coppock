// -------------------------------------------------------------------------------
//   Classical Coppock indicator with multi-timeframe support.
//   Should be applied to monthly timeframe for classic usage.
//   Change timeframe and parameters for experimental usage.
//   
//   Version 1.01
//   Copyright 2025, EarnForex.com
//   https://www.earnforex.com/indicators/Coppock/
// -------------------------------------------------------------------------------

using System;
using System.Collections.Generic;
using cAlgo.API;
using cAlgo.API.Indicators;
using cAlgo.API.Internals;

namespace cAlgo.Indicators
{
    [Indicator(IsOverlay = false, TimeZone = TimeZones.UTC, AccessRights = AccessRights.None)]
    public class Coppock : Indicator
    {
        public enum CandleToCheck
        {
            CurrentCandle,
            ClosedCandle
        }

        public enum TradeSignal
        {
            Buy = 1,
            Sell = -1,
            Neutral = 0
        }

        [Parameter("Timeframe", DefaultValue = "Current")]
        public TimeFrame InputTimeFrame { get; set; }

        [Parameter("ROC1 Period", DefaultValue = 14, MinValue = 1)]
        public int ROC1Period { get; set; }

        [Parameter("ROC2 Period", DefaultValue = 11, MinValue = 1)]
        public int ROC2Period { get; set; }

        [Parameter("MA Period", DefaultValue = 10, MinValue = 1)]
        public int MAPeriod { get; set; }

        [Parameter("MA Type", DefaultValue = MovingAverageType.Weighted)]
        public MovingAverageType MAType { get; set; }

        [Parameter("Strict Signals", DefaultValue = true)]
        public bool StrictSignals { get; set; }

        [Parameter("Trigger Candle", DefaultValue = CandleToCheck.CurrentCandle)]
        public CandleToCheck TriggerCandle { get; set; }

        [Parameter("Enable Native Alerts", DefaultValue = false)]
        public bool EnableNativeAlerts { get; set; }

        [Parameter("Enable Email Alerts", DefaultValue = false)]
        public bool EnableEmailAlerts { get; set; }

        [Parameter("Email Address", DefaultValue = "", Group = "Alerts")]
        public string EmailAddress { get; set; }

        [Parameter("Enable Sound Alerts", DefaultValue = false)]
        public bool EnableSoundAlerts { get; set; }

        [Parameter("Sound Type", DefaultValue = SoundType.Announcement, Group = "Alerts")]
        public SoundType SoundType { get; set; }

        [Parameter("Wait Time Between Notifications (Seconds)", DefaultValue = 5, MinValue = 1)]
        public int WaitTimeNotify { get; set; }

        [Parameter("Draw Signal Arrows", DefaultValue = true)]
        public bool EnableDrawArrows { get; set; }

        [Parameter("Arrow Size", DefaultValue = 3, MinValue = 1, MaxValue = 5)]
        public int ArrowSize { get; set; }

        [Parameter("Buy Arrow Color", DefaultValue = "Green")]
        public Color ArrowBuyColor { get; set; }

        [Parameter("Sell Arrow Color", DefaultValue = "Red")]
        public Color ArrowSellColor { get; set; }

        [Parameter("Indicator Short Name", DefaultValue = "COP")]
        public string IndicatorName { get; set; }

        [Output("Coppock", LineColor = "Red", LineStyle = LineStyle.Solid, PlotType = PlotType.Histogram, Thickness = 2)]
        public IndicatorDataSeries CoppockOutput { get; set; }

        private IndicatorDataSeries ROCSum;
        private MovingAverage MA;
        private Bars MTFBars;
        private DateTime LastNotificationTime = DateTime.MinValue;
        private TradeSignal LastSignal = TradeSignal.Neutral;
        private int DrawBegin;
        private bool IsMTF;
        private Dictionary<int, double> maCache = new Dictionary<int, double>();
        private Dictionary<int, double> ma2Cache = new Dictionary<int, double>();        
        private Dictionary<int, double> ma3Cache = new Dictionary<int, double>();        

        protected override void Initialize()
        {
            // Initialize MTF bars.
            if (InputTimeFrame != null && InputTimeFrame != Chart.TimeFrame)
            {
                MTFBars = MarketData.GetBars(InputTimeFrame);
                IsMTF = true;
            }
            else
            {
                MTFBars = Bars;
                IsMTF = false;
            }

            // Calculate draw begin.
            DrawBegin = Math.Max(ROC1Period, ROC2Period) + MAPeriod;

            // Initialize ROCSum buffer.
            ROCSum = CreateDataSeries();

            if (!IsMTF)
            {
                // MA will be applied directly to ROCSum for non-MTF.
                MA = Indicators.MovingAverage(ROCSum, MAPeriod, MAType);
            }

            // Set indicator name.
            string maTypeStr = MAType.ToString();
            string shortName = $"Coppock({ROC1Period}, {ROC2Period}) {maTypeStr}({MAPeriod})";
            if (MTFBars != Bars)
            {
                shortName += $" {InputTimeFrame}";
            }
        }

        public override void Calculate(int index)
        {
            // Skip if not enough bars.
            if (index < DrawBegin)
            {
                CoppockOutput[index] = 0;
                return;
            }

            // Calculate based on selected timeframe.
            if (IsMTF)
            {
                CalculateMTF(index);
            }
            else
            {
                CalculateNormal(index);
            }

            // Check for signals and draw arrows on the last bar.
            if (EnableDrawArrows)
            {
                if (!IsMTF || !IsLastBar) CheckAndDrawArrow(index);
                else
                {
                    int currentMTFIndex = GetMTFIndex(index);
                    if (currentMTFIndex < 0) return; // Couldn't find MTF index.
                    for (int j = index - 1; j >= DrawBegin; j--)
                    {
                        int prevMTFIndex = GetMTFIndex(j);
                        if (prevMTFIndex < 0) return; // Couldn't find MTF index.
                        if (prevMTFIndex != currentMTFIndex)
                        {
                            CheckAndDrawArrow(j + 1); // Run only on the first LTF bar of the current HTF bar.
                            break;
                        }
                    }
                }
            }

            // Check alerts.
            if (IsLastBar && (EnableNativeAlerts || EnableEmailAlerts || EnableSoundAlerts))
            {
                if (!IsMTF) CheckAlerts(index);
                else
                {
                    int currentMTFIndex = GetMTFIndex(index);
                    if (currentMTFIndex < 0) return; // Couldn't find MTF index.
                    for (int j = index - 1; j >= DrawBegin; j--)
                    {
                        int prevMTFIndex = GetMTFIndex(j);
                        if (prevMTFIndex < 0) return; // Couldn't find MTF index.
                        if (prevMTFIndex != currentMTFIndex)
                        {
                            CheckAlerts(j + 1); // Run only on the first LTF bar of the current HTF bar.
                            break;
                        }
                    }
                }
            }
        }

        private void CalculateNormal(int index)
        {
            // Calculate ROC1 and ROC2.
            if (index >= ROC1Period && index >= ROC2Period)
            {
                double currentClose = Bars.ClosePrices[index];
                double closeROC1 = Bars.ClosePrices[index - ROC1Period];
                double closeROC2 = Bars.ClosePrices[index - ROC2Period];

                if (closeROC1 != 0 && closeROC2 != 0)
                {
                    double roc1 = (currentClose - closeROC1) / closeROC1;
                    double roc2 = (currentClose - closeROC2) / closeROC2;
                    ROCSum[index] = roc1 + roc2;
                }
                else
                {
                    ROCSum[index] = 0;
                }
            }
            else
            {
                ROCSum[index] = 0;
            }

            // Apply moving average.
            CoppockOutput[index] = MA.Result[index];
        }

        private void CalculateMTF(int index)
        {
            // Find corresponding MTF bar index.
            int mtfIndex = GetMTFIndex(index);
            
            if (mtfIndex < 0 || mtfIndex >= MTFBars.Count)
            {
                CoppockOutput[index] = 0;
                return;
            }

            int mtfIndex_prev = mtfIndex;
            while (mtfIndex_prev == mtfIndex) // Need to check (or re-check) all LTF bars of the current HTF bar.
            {
                // Calculate ROC values and MA for this MTF bar if not already done.
                if (mtfIndex >= Math.Max(ROC1Period, ROC2Period))
                {
                    double currentClose = MTFBars.ClosePrices[mtfIndex];
                    double closeROC1 = MTFBars.ClosePrices[mtfIndex - ROC1Period];
                    double closeROC2 = MTFBars.ClosePrices[mtfIndex - ROC2Period];
                    
                    double rocSum = 0;
                    if (closeROC1 != 0 && closeROC2 != 0)
                    {
                        rocSum = (currentClose - closeROC1) / closeROC1 + 
                                (currentClose - closeROC2) / closeROC2;
                    }
                    
                    // Calculate MA manually for MTF data.
                    if (mtfIndex >= DrawBegin)
                    {
                        double maValue = CalculateMAForMTF(mtfIndex);
                        CoppockOutput[index] = maValue;
                    }
                    else
                    {
                        CoppockOutput[index] = 0;
                    }
                }
                else
                {
                    CoppockOutput[index] = 0;
                }
                if (IsLastBar) // This has to be done only for the current (latest) bar.
                {
                    index--; // Go back one LTF bar.
                    mtfIndex = GetMTFIndex(index); // Retrieve new HTF bar.
                }
                else
                {
                    break; // If this is one of the historical bars, no need to cycle back.
                }
            }
        }

        private double CalculateMAForMTF(int mtfIndex)
        {
            // Calculate ROCSum values for the MA window.
            double[] rocValues = new double[MAPeriod];
            
            for (int i = 0; i < MAPeriod; i++)
            {
                int idx = mtfIndex - i;
                if (idx >= Math.Max(ROC1Period, ROC2Period))
                {
                    double currentClose = MTFBars.ClosePrices[idx];
                    double closeROC1 = MTFBars.ClosePrices[idx - ROC1Period];
                    double closeROC2 = MTFBars.ClosePrices[idx - ROC2Period];
                    
                    if (closeROC1 != 0 && closeROC2 != 0)
                    {
                        rocValues[i] = (currentClose - closeROC1) / closeROC1 + 
                                      (currentClose - closeROC2) / closeROC2;
                    }
                }
            }
            
            // Calculate the appropriate MA type.
            switch (MAType)
            {
                case MovingAverageType.Simple:
                    return CalculateSMA(rocValues);
                    
                case MovingAverageType.Exponential:
                    return CalculateEMA(mtfIndex);
                    
                case MovingAverageType.Weighted:
                    return CalculateLWMA(rocValues);
                    
                case MovingAverageType.WilderSmoothing:
                    return CalculateSMMA(mtfIndex);
                
                case MovingAverageType.Hull:
                    return CalculateHullMA(mtfIndex);
                
                case MovingAverageType.Triangular:
                    return CalculateTriangularMA(mtfIndex);
                
                case MovingAverageType.VIDYA:
                    return CalculateVIDYA(mtfIndex);

                case MovingAverageType.TimeSeries:
                    return CalculateTimeSeriesMA(mtfIndex);

                case MovingAverageType.KaufmanAdaptive:
                    return CalculateKAMA(mtfIndex);
                    
                case MovingAverageType.DoubleExponential:
                    return CalculateDEMA(mtfIndex);

                case MovingAverageType.TripleExponential:
                    return CalculateTEMA(mtfIndex);

                default:
                    return CalculateSMA(rocValues);
            }
        }
        
        private double CalculateSMA(double[] values)
        {
            double sum = 0;
            foreach (var value in values)
                sum += value;
            return sum / values.Length;
        }
        
        private double CalculateLWMA(double[] values)
        {
            double sum = 0;
            double weightSum = 0;
            for (int i = 0; i < values.Length; i++)
            {
                int weight = values.Length - i;
                sum += values[i] * weight;
                weightSum += weight;
            }
            return sum / weightSum;
        }
        
        private double CalculateEMA(int mtfIndex)
        {
            double alpha = 2.0 / (MAPeriod + 1.0);
            
            // If we have a cached value for the previous bar, use it.
            if (maCache.ContainsKey(mtfIndex - 1))
            {
                double currentROC = 0;
                if (mtfIndex >= Math.Max(ROC1Period, ROC2Period))
                {
                    double currentClose = MTFBars.ClosePrices[mtfIndex];
                    double closeROC1 = MTFBars.ClosePrices[mtfIndex - ROC1Period];
                    double closeROC2 = MTFBars.ClosePrices[mtfIndex - ROC2Period];
                    
                    if (closeROC1 != 0 && closeROC2 != 0)
                    {
                        currentROC = (currentClose - closeROC1) / closeROC1 + 
                                     (currentClose - closeROC2) / closeROC2;
                    }
                }
                
                double ema = alpha * currentROC + (1 - alpha) * maCache[mtfIndex - 1];
                maCache[mtfIndex] = ema;
                return ema;
            }
            else
            {
                // Initialize with SMA for the first calculation.
                double[] rocValues = new double[MAPeriod];
                for (int i = 0; i < MAPeriod; i++)
                {
                    int idx = mtfIndex - i;
                    if (idx >= Math.Max(ROC1Period, ROC2Period))
                    {
                        double currentClose = MTFBars.ClosePrices[idx];
                        double closeROC1 = MTFBars.ClosePrices[idx - ROC1Period];
                        double closeROC2 = MTFBars.ClosePrices[idx - ROC2Period];
                        
                        if (closeROC1 != 0 && closeROC2 != 0)
                        {
                            rocValues[i] = (currentClose - closeROC1) / closeROC1 + 
                                           (currentClose - closeROC2) / closeROC2;
                        }
                    }
                }
                double sma = CalculateSMA(rocValues);
                maCache[mtfIndex] = sma;
                return sma;
            }
        }

        private double CalculateSMMA(int mtfIndex)
        {
            // Wilder's Smoothing uses alpha = 1/period instead of 2/(period+1) like EMA.
            double alpha = 1.0 / MAPeriod;
            
            // If we have a cached value for the previous bar, use it.
            if (maCache.ContainsKey(mtfIndex - 1))
            {
                double currentROC = 0;
                if (mtfIndex >= Math.Max(ROC1Period, ROC2Period))
                {
                    double currentClose = MTFBars.ClosePrices[mtfIndex];
                    double closeROC1 = MTFBars.ClosePrices[mtfIndex - ROC1Period];
                    double closeROC2 = MTFBars.ClosePrices[mtfIndex - ROC2Period];
                    
                    if (closeROC1 != 0 && closeROC2 != 0)
                    {
                        currentROC = (currentClose - closeROC1) / closeROC1 + 
                                     (currentClose - closeROC2) / closeROC2;
                    }
                }
                
                // Wilder's formula: WMA[i] = (Previous WMA * (n-1) + Data[i]) / n
                // Which is equivalent to: WMA[i] = alpha * Data[i] + (1 - alpha) * Previous WMA.
                double wilder = alpha * currentROC + (1 - alpha) * maCache[mtfIndex - 1];
                maCache[mtfIndex] = wilder;
                return wilder;
            }
            else
            {
                // Initialize with SMA for the first calculation (standard practice for Wilder's).
                double[] rocValues = new double[MAPeriod];
                for (int i = 0; i < MAPeriod; i++)
                {
                    int idx = mtfIndex - i;
                    if (idx >= Math.Max(ROC1Period, ROC2Period))
                    {
                        double currentClose = MTFBars.ClosePrices[idx];
                        double closeROC1 = MTFBars.ClosePrices[idx - ROC1Period];
                        double closeROC2 = MTFBars.ClosePrices[idx - ROC2Period];
                        
                        if (closeROC1 != 0 && closeROC2 != 0)
                        {
                            rocValues[i] = (currentClose - closeROC1) / closeROC1 + 
                                           (currentClose - closeROC2) / closeROC2;
                        }
                    }
                }
                double sma = CalculateSMA(rocValues);
                maCache[mtfIndex] = sma;
                return sma;
            }
        }

        private double CalculateHullMA(int mtfIndex)
        {
            // Hull Moving Average formula:
            // HMA = WMA(2*WMA(n/2) - WMA(n), sqrt(n))
            // Where n is the period.
            
            int halfPeriod = MAPeriod / 2;
            int sqrtPeriod = (int)Math.Sqrt(MAPeriod);
            
            // Need enough bars for the calculation.
            if (mtfIndex < MAPeriod + sqrtPeriod - 1)
                return 0;
            
            // Step 1: Calculate WMA with period n/2.
            double wmaHalf = CalculateWMA_MTF(mtfIndex, halfPeriod);
            
            // Step 2: Calculate WMA with full period n.
            double wmaFull = CalculateWMA_MTF(mtfIndex, MAPeriod);
            
            // Step 3: Calculate the difference series: 2*WMA(n/2) - WMA(n).
            // We need to calculate this for sqrt(n) bars.
            double[] diffSeries = new double[sqrtPeriod];
            for (int i = 0; i < sqrtPeriod; i++)
            {
                int idx = mtfIndex - i;
                if (idx >= MAPeriod - 1)
                {
                    double wmaHalfAtIdx = CalculateWMA_MTF(idx, halfPeriod);
                    double wmaFullAtIdx = CalculateWMA_MTF(idx, MAPeriod);
                    diffSeries[i] = 2 * wmaHalfAtIdx - wmaFullAtIdx;
                }
                else
                {
                    diffSeries[i] = 0;
                }
            }
            
            // Step 4: Calculate final WMA with period sqrt(n) on the difference series.
            double sum = 0;
            double weightSum = 0;
            for (int i = 0; i < sqrtPeriod; i++)
            {
                int weight = sqrtPeriod - i;
                sum += diffSeries[i] * weight;
                weightSum += weight;
            }
            
            return weightSum > 0 ? sum / weightSum : 0;
        }

        // Helper method to calculate Weighted Moving Average for a specific period.
        private double CalculateWMA_MTF(int mtfIndex, int period)
        {
            if (mtfIndex < period - 1) return 0;
                
            double sum = 0;
            double weightSum = 0;
            
            for (int i = 0; i < period; i++)
            {
                int idx = mtfIndex - i;
                if (idx >= Math.Max(ROC1Period, ROC2Period))
                {
                    double currentClose = MTFBars.ClosePrices[idx];
                    double closeROC1 = MTFBars.ClosePrices[idx - ROC1Period];
                    double closeROC2 = MTFBars.ClosePrices[idx - ROC2Period];
                    
                    if (closeROC1 != 0 && closeROC2 != 0)
                    {
                        double rocValue = (currentClose - closeROC1) / closeROC1 + 
                                          (currentClose - closeROC2) / closeROC2;
                        int weight = period - i;
                        sum += rocValue * weight;
                        weightSum += weight;
                    }
                }
            }
            
            return weightSum > 0 ? sum / weightSum : 0;
        }

        private double CalculateTriangularMA(int mtfIndex)
        {
            // Triangular Moving Average is a double-smoothed SMA.
            // TMA = SMA(SMA(data, ceil(n/2)), floor(n/2)+1)
            // Or alternatively: TMA = SMA(SMA(data, period1), period2)
            // Where period1 = ceil(period/2) and period2 = floor(period/2) + 1.
            
            int period1 = (MAPeriod + 1) / 2;  // Ceiling of MAPeriod/2.
            int period2 = MAPeriod / 2 + 1;    // Floor of MAPeriod/2 + 1.
            
            // Need enough bars for the calculation.
            if (mtfIndex < period1 + period2 - 2)
                return 0;
            
            // Step 1: Calculate first SMA with period1.
            // We need period2 values of this first SMA to calculate the final TMA.
            double[] firstSMA = new double[period2];
            
            for (int i = 0; i < period2; i++)
            {
                int idx = mtfIndex - i;
                if (idx >= period1 - 1)
                {
                    firstSMA[i] = CalculateSMA_MTF(idx, period1);
                }
                else
                {
                    firstSMA[i] = 0;
                }
            }
            
            // Step 2: Calculate second SMA on the first SMA values.
            double sum = 0;
            int validCount = 0;
            for (int i = 0; i < period2; i++)
            {
                if (firstSMA[i] != 0)
                {
                    sum += firstSMA[i];
                    validCount++;
                }
            }
            
            return validCount > 0 ? sum / validCount : 0;
        }
        
        private double CalculateSMA_MTF(int mtfIndex, int period)
        {
            // Helper method to calculate Simple Moving Average for a specific period on MTF data.
            if (mtfIndex < period - 1) return 0;
                
            double sum = 0;
            int validCount = 0;
            
            for (int i = 0; i < period; i++)
            {
                int idx = mtfIndex - i;
                if (idx >= Math.Max(ROC1Period, ROC2Period))
                {
                    double currentClose = MTFBars.ClosePrices[idx];
                    double closeROC1 = MTFBars.ClosePrices[idx - ROC1Period];
                    double closeROC2 = MTFBars.ClosePrices[idx - ROC2Period];
                    
                    if (closeROC1 != 0 && closeROC2 != 0)
                    {
                        double rocValue = (currentClose - closeROC1) / closeROC1 + 
                                          (currentClose - closeROC2) / closeROC2;
                        sum += rocValue;
                        validCount++;
                    }
                }
            }
            
            return validCount > 0 ? sum / validCount : 0;
        }

        private double CalculateVIDYA(int mtfIndex)
        {
            // VIDYA (Variable Index Dynamic Average) uses Chande Momentum Oscillator (CMO) for volatility.
            // Formula: VIDYA = α × VI × Close + (1 - α × VI) × Previous VIDYA
            // Where:
            // α = 2 / (Period + 1) - standard EMA smoothing constant.
            // VI = |CMO| / 100 - Volatility Index based on Chande Momentum Oscillator.
            // CMO = (Sum of positive changes - Sum of negative changes) / (Sum of positive changes + Sum of negative changes) × 100.
            
            // We'll use a 9-period CMO as standard (can be parameterized if needed).
            int cmoPeriod = 9;
            double alpha = 2.0 / (MAPeriod + 1.0);
            
            // Need enough bars for CMO calculation.
            if (mtfIndex < Math.Max(cmoPeriod, MAPeriod))
                return 0;
            
            // Calculate CMO (Chande Momentum Oscillator) for volatility index.
            double sumUp = 0;
            double sumDown = 0;
            
            for (int i = 0; i < cmoPeriod; i++)
            {
                int currentIdx = mtfIndex - i;
                int prevIdx = currentIdx - 1;
                
                if (currentIdx >= Math.Max(ROC1Period, ROC2Period) && 
                    prevIdx >= Math.Max(ROC1Period, ROC2Period))
                {
                    // Calculate ROC values for current and previous bars.
                    double currROC = CalculateROC_MTF(currentIdx);
                    double prevROC = CalculateROC_MTF(prevIdx);
                    
                    double change = currROC - prevROC;
                    if (change > 0)
                        sumUp += change;
                    else
                        sumDown += Math.Abs(change);
                }
            }
            
            // Calculate CMO and Volatility Index.
            double cmo = 0;
            if (sumUp + sumDown != 0)
            {
                cmo = ((sumUp - sumDown) / (sumUp + sumDown)) * 100;
            }
            double vi = Math.Abs(cmo) / 100.0; // Volatility Index (0 to 1).
            
            // Get current ROC value.
            double currentROC = CalculateROC_MTF(mtfIndex);
            
            // Calculate VIDYA.
            if (maCache.ContainsKey(mtfIndex - 1))
            {
                // VIDYA formula with dynamic smoothing.
                double vidya = (alpha * vi * currentROC) + (1 - alpha * vi) * maCache[mtfIndex - 1];
                maCache[mtfIndex] = vidya;
                return vidya;
            }
            else
            {
                // Initialize with SMA for the first calculation.
                double sma = CalculateSMA_MTF(mtfIndex, MAPeriod);
                maCache[mtfIndex] = sma;
                return sma;
            }
        }
        
        private double CalculateROC_MTF(int mtfIndex)
        {
            // Helper method to calculate ROC value for a specific MTF bar.
            if (mtfIndex < Math.Max(ROC1Period, ROC2Period))
                return 0;
                
            double currentClose = MTFBars.ClosePrices[mtfIndex];
            double closeROC1 = MTFBars.ClosePrices[mtfIndex - ROC1Period];
            double closeROC2 = MTFBars.ClosePrices[mtfIndex - ROC2Period];
            
            if (closeROC1 != 0 && closeROC2 != 0)
            {
                return (currentClose - closeROC1) / closeROC1 + 
                       (currentClose - closeROC2) / closeROC2;
            }
            
            return 0;
        }

        private double CalculateTimeSeriesMA(int mtfIndex)
        {
            // Time Series Moving Average is based on linear regression.
            // It calculates the linear regression line over the period and returns the current predicted value.
            // Formula: TSMA = a + b * n
            // Where:
            // a = intercept of the regression line.
            // b = slope of the regression line.
            // n = period (we use n = period for the endpoint of the regression line).
            // 
            // Linear regression formulas:
            // b = (n * Σ(x*y) - Σx * Σy) / (n * Σ(x²) - (Σx)²)
            // a = (Σy - b * Σx) / n
            
            // Need enough bars for the calculation.
            if (mtfIndex < MAPeriod - 1)
                return 0;
            
            // Prepare data for linear regression.
            double sumX = 0;      // Sum of x values (indices).
            double sumY = 0;      // Sum of y values (ROC values).
            double sumXY = 0;     // Sum of x*y.
            double sumX2 = 0;     // Sum of x².
            int validCount = 0;
            
            for (int i = 0; i < MAPeriod; i++)
            {
                int idx = mtfIndex - (MAPeriod - 1 - i);  // Start from oldest bar in the period.
                if (idx >= Math.Max(ROC1Period, ROC2Period))
                {
                    double rocValue = CalculateROC_MTF(idx);
                    
                    // x values are 1, 2, 3, ..., MAPeriod.
                    double x = i + 1;
                    double y = rocValue;
                    
                    sumX += x;
                    sumY += y;
                    sumXY += x * y;
                    sumX2 += x * x;
                    validCount++;
                }
            }
            
            if (validCount < 2)  // Need at least 2 points for regression.
                return 0;
            
            // Calculate slope (b) and intercept (a).
            double denominator = validCount * sumX2 - sumX * sumX;
            if (Math.Abs(denominator) < 0.0000001)  // Avoid division by zero.
                return validCount > 0 ? sumY / validCount : 0;  // Return simple average if no slope.
            
            double slope = (validCount * sumXY - sumX * sumY) / denominator;
            double intercept = (sumY - slope * sumX) / validCount;
            
            // Time Series MA returns the predicted value at the end of the regression line.
            // This is the value at x = MAPeriod.
            double tsma = intercept + slope * MAPeriod;
            
            return tsma;
        }

        private double CalculateKAMA(int mtfIndex)
        {
            // Kaufman Adaptive Moving Average (KAMA).
            // Adjusts its sensitivity based on the Efficiency Ratio (ER).
            // Formula: KAMA = Previous KAMA + SC × (Price - Previous KAMA)
            // Where SC = [ER × (FastSC - SlowSC) + SlowSC]²
            // ER = Change / Volatility.
            
            // KAMA standard parameters.
            int erPeriod = MAPeriod;  // Efficiency Ratio period (typically 10).
            int fastPeriod = 2;       // Fast EMA period (typically 2).
            int slowPeriod = 30;      // Slow EMA period (typically 30).
            
            // Need enough bars for ER calculation.
            if (mtfIndex < erPeriod) return 0;
            
            // Calculate Efficiency Ratio (ER).
            // ER = Direction / Volatility.
            // Direction = ABS(Close - Close[n]).
            // Volatility = Sum of ABS(Close - Close[1]) over n periods.
            
            double direction = 0;
            double volatility = 0;
            
            // Get current ROC value.
            double currentROC = CalculateROC_MTF(mtfIndex);
            
            if (mtfIndex >= erPeriod + Math.Max(ROC1Period, ROC2Period))
            {
                // Direction: absolute change over the period.
                double oldROC = CalculateROC_MTF(mtfIndex - erPeriod);
                direction = Math.Abs(currentROC - oldROC);
                
                // Volatility: sum of absolute changes.
                for (int i = 0; i < erPeriod; i++)
                {
                    double roc1 = CalculateROC_MTF(mtfIndex - i);
                    double roc2 = CalculateROC_MTF(mtfIndex - i - 1);
                    volatility += Math.Abs(roc1 - roc2);
                }
            }
            
            // Calculate Efficiency Ratio.
            double er = 0;
            if (volatility != 0)
            {
                er = direction / volatility;
            }
            
            // Calculate Smoothing Constant (SC).
            // Fast SC = 2 / (fastPeriod + 1).
            // Slow SC = 2 / (slowPeriod + 1).
            double fastSC = 2.0 / (fastPeriod + 1.0);
            double slowSC = 2.0 / (slowPeriod + 1.0);
            
            // SC = [ER × (fastSC - slowSC) + slowSC]².
            double sc = er * (fastSC - slowSC) + slowSC;
            sc = sc * sc;  // Square the result.
            
            // Calculate KAMA.
            if (maCache.ContainsKey(mtfIndex - 1))
            {
                // KAMA formula.
                double kama = maCache[mtfIndex - 1] + sc * (currentROC - maCache[mtfIndex - 1]);
                maCache[mtfIndex] = kama;
                return kama;
            }
            else
            {
                // Initialize with SMA for the first calculation.
                double sma = CalculateSMA_MTF(mtfIndex, MAPeriod);
                maCache[mtfIndex] = sma;
                return sma;
            }
        }

        private double CalculateDEMA(int mtfIndex)
        {
            // Double Exponential Moving Average (DEMA).
            // DEMA = 2 * EMA - EMA(EMA)
            // This reduces the lag of traditional EMA by compensating for the smoothing.
            
            // Need enough bars for the calculation.
            if (mtfIndex < MAPeriod - 1) return 0;
            
            double alpha = 2.0 / (MAPeriod + 1.0);
            
            // Step 1: Calculate first EMA of ROC values.
            double ema1;
            double currentROC = CalculateROC_MTF(mtfIndex);
            
            if (maCache.ContainsKey(mtfIndex - 1))
            {
                // Continue with EMA calculation.
                ema1 = alpha * currentROC + (1 - alpha) * maCache[mtfIndex - 1];
            }
            else
            {
                // Initialize with SMA.
                double sum = 0;
                int count = 0;
                for (int i = 0; i < MAPeriod; i++)
                {
                    int idx = mtfIndex - i;
                    if (idx >= Math.Max(ROC1Period, ROC2Period))
                    {
                        sum += CalculateROC_MTF(idx);
                        count++;
                    }
                }
                ema1 = count > 0 ? sum / count : 0;
            }
            maCache[mtfIndex] = ema1;
            
            // Step 2: Calculate EMA of EMA (second smoothing).
            double ema2;
            
            if (ma2Cache.ContainsKey(mtfIndex - 1))
            {
                // Use the first EMA as input for the second EMA.
                ema2 = alpha * ema1 + (1 - alpha) * ma2Cache[mtfIndex - 1];
            }
            else
            {
                // Initialize second EMA with SMA of first EMA values.
                double sum = 0;
                int count = 0;
                
                // We need to calculate historical EMA1 values if not cached.
                for (int i = MAPeriod - 1; i >= 0; i--)
                {
                    int idx = mtfIndex - i;
                    if (idx >= MAPeriod - 1)
                    {
                        // Try to get from cache or calculate.
                        if (maCache.ContainsKey(idx))
                        {
                            sum += maCache[idx];
                            count++;
                        }
                        else if (idx == mtfIndex)
                        {
                            sum += ema1;
                            count++;
                        }
                    }
                }
                
                ema2 = count > 0 ? sum / count : ema1;
            }
            ma2Cache[mtfIndex] = ema2;
            
            // Step 3: Calculate DEMA = 2 * EMA - EMA(EMA).
            double dema = 2 * ema1 - ema2;
            
            return dema;
        }

        private double CalculateTEMA(int mtfIndex)
        {
            // Triple Exponential Moving Average (TEMA).
            // TEMA = 3 * EMA - 3 * EMA(EMA) + EMA(EMA(EMA))
            // This further reduces lag compared to DEMA by adding a third smoothing level.
            
            // Need enough bars for the calculation.
            if (mtfIndex < MAPeriod - 1)
                return 0;
            
            double alpha = 2.0 / (MAPeriod + 1.0);
            
            // Step 1: Calculate first EMA of ROC values.
            double ema1;
            double currentROC = CalculateROC_MTF(mtfIndex);
            
            if (maCache.ContainsKey(mtfIndex - 1))
            {
                // Continue with EMA calculation.
                ema1 = alpha * currentROC + (1 - alpha) * maCache[mtfIndex - 1];
            }
            else
            {
                // Initialize with SMA.
                double sum = 0;
                int count = 0;
                for (int i = 0; i < MAPeriod; i++)
                {
                    int idx = mtfIndex - i;
                    if (idx >= Math.Max(ROC1Period, ROC2Period))
                    {
                        sum += CalculateROC_MTF(idx);
                        count++;
                    }
                }
                ema1 = count > 0 ? sum / count : 0;
            }
            maCache[mtfIndex] = ema1;
            
            // Step 2: Calculate EMA of EMA (second smoothing).
            double ema2;
            
            if (ma2Cache.ContainsKey(mtfIndex - 1))
            {
                // Use the first EMA as input for the second EMA.
                ema2 = alpha * ema1 + (1 - alpha) * ma2Cache[mtfIndex - 1];
            }
            else
            {
                // Initialize second EMA with SMA of first EMA values.
                double sum = 0;
                int count = 0;
                
                // Calculate historical EMA1 values if needed.
                for (int i = MAPeriod - 1; i >= 0; i--)
                {
                    int idx = mtfIndex - i;
                    if (idx >= MAPeriod - 1)
                    {
                        if (maCache.ContainsKey(idx))
                        {
                            sum += maCache[idx];
                            count++;
                        }
                        else if (idx == mtfIndex)
                        {
                            sum += ema1;
                            count++;
                        }
                    }
                }
                
                ema2 = count > 0 ? sum / count : ema1;
            }
            ma2Cache[mtfIndex] = ema2;
            
            // Step 3: Calculate EMA of EMA of EMA (third smoothing).
            double ema3;
            
            if (ma3Cache.ContainsKey(mtfIndex - 1))
            {
                // Use the second EMA as input for the third EMA.
                ema3 = alpha * ema2 + (1 - alpha) * ma3Cache[mtfIndex - 1];
            }
            else
            {
                // Initialize third EMA with SMA of second EMA values.
                double sum = 0;
                int count = 0;
                
                // Calculate historical EMA2 values if needed.
                for (int i = MAPeriod - 1; i >= 0; i--)
                {
                    int idx = mtfIndex - i;
                    if (idx >= MAPeriod - 1)
                    {
                        if (ma2Cache.ContainsKey(idx))
                        {
                            sum += ma2Cache[idx];
                            count++;
                        }
                        else if (idx == mtfIndex)
                        {
                            sum += ema2;
                            count++;
                        }
                    }
                }
                
                ema3 = count > 0 ? sum / count : ema2;
            }
            ma3Cache[mtfIndex] = ema3;
            
            // Step 4: Calculate TEMA = 3 * EMA - 3 * EMA(EMA) + EMA(EMA(EMA)).
            double tema = 3 * ema1 - 3 * ema2 + ema3;
            
            return tema;
        }

        private int GetMTFIndex(int currentIndex)
        {
            if (MTFBars == Bars)
                return currentIndex;

            DateTime currentTime = Bars.OpenTimes[currentIndex];
            
            // Binary search for efficiency.
            int left = 0;
            int right = MTFBars.Count - 1;
            
            while (left <= right)
            {
                int mid = (left + right) / 2;
                DateTime mtfTime = MTFBars.OpenTimes[mid];
                
                if (currentTime >= mtfTime && 
                    (mid == MTFBars.Count - 1 || currentTime < MTFBars.OpenTimes[mid + 1]))
                {
                    return mid;
                }
                
                if (currentTime < mtfTime)
                {
                    right = mid - 1;
                }
                else
                {
                    left = mid + 1;
                }
            }
            
            return -1;
        }

        private TradeSignal IsSignal(int index)
        {
            int shift = TriggerCandle == CandleToCheck.ClosedCandle ? 1 : 0;
            int checkIndex = index;
            
            if (!IsMTF) checkIndex -= shift;
            else if (TriggerCandle == CandleToCheck.ClosedCandle) // Find the first LTF candle inside the previous MTF index when working with Previous Candle signals.
            {
                int currentMTFIndex = GetMTFIndex(index);
                if (currentMTFIndex < 0) return TradeSignal.Neutral; // Couldn't find MTF index.
                for (int j = index - shift; j >= DrawBegin; j--)
                {
                    int prevMTFIndex = GetMTFIndex(j);
                    if (prevMTFIndex < 0) return TradeSignal.Neutral; // Couldn't find MTF index.
                    if (prevMTFIndex != currentMTFIndex)
                    {
                        for (int k = j - 1; k >= DrawBegin; k--)
                        {
                            int prevprevMTFIndex = GetMTFIndex(k);
                            if (prevprevMTFIndex < 0) return TradeSignal.Neutral; // Couldn't find MTF index.
                            if (prevprevMTFIndex != prevMTFIndex)
                            {
                                checkIndex = k + 1; // Returning to the previous step because overstepped to the earlier HTF bar.
                                break;
                            }
                        }
                        break;
                    }
                }
            }

            if (checkIndex <= DrawBegin + 1)
                return TradeSignal.Neutral;

            // Find previous direction.
            int prevDirection = 0;
            for (int i = checkIndex - 2; i >= DrawBegin; i--)
            {
                if (CoppockOutput[i] > CoppockOutput[checkIndex - 1])
                {
                    prevDirection = -1; // Was falling.
                    break;
                }
                else if (CoppockOutput[i] < CoppockOutput[checkIndex - 1])
                {
                    prevDirection = 1; // Was rising.
                    break;
                }
            }

            if (prevDirection == 0)
                return TradeSignal.Neutral;

            // Check current direction change.
            if (CoppockOutput[checkIndex] > CoppockOutput[checkIndex - 1])
            {
                // Now rising.
                if (prevDirection == -1) // Was falling, now rising = Buy signal.
                {
                    if (!StrictSignals || CoppockOutput[checkIndex] < 0)
                        return TradeSignal.Buy;
                }
            }
            else if (CoppockOutput[checkIndex] < CoppockOutput[checkIndex - 1])
            {
                // Now falling.
                if (prevDirection == 1) // Was rising, now falling = Sell signal.
                {
                    if (!StrictSignals || CoppockOutput[checkIndex] > 0)
                        return TradeSignal.Sell;
                }
            }

            return TradeSignal.Neutral;
        }

        private void CheckAndDrawArrow(int index)
        {
            TradeSignal signal = IsSignal(index);
            
            if (signal == TradeSignal.Neutral)
            {
                if (IsLastBar || (IsMTF && IsLastBar && GetMTFIndex(index) == MTFBars.Count - 1)) RemoveCurrentArrow(index); // Call up only for the latest bar (bars in MTF).
            }
            else
            {
                if (IsMTF && (index == 0 || GetMTFIndex(index) == GetMTFIndex(index - 1))) return; // Avoid putting an arrow on all the LTF bars of a trigger HTF candle in MTF mode.
                DrawArrowObject(signal, index);
            }
        }

        private void DrawArrowObject(TradeSignal signal, int index)
        {
            DateTime arrowTime = Bars.OpenTimes[index];
            string arrowName = $"{IndicatorName}-ARROW-{arrowTime:yyyyMMddHHmmss}";

            // Remove existing arrow at this position.
            Chart.RemoveObject(arrowName);

            double arrowPrice;
            string arrowText;
            Color arrowColor;
            VerticalAlignment vAlign;

            if (signal == TradeSignal.Buy)
            {
                arrowPrice = Bars.LowPrices[index];
                arrowText = "▲";
                arrowColor = ArrowBuyColor;
                vAlign = VerticalAlignment.Bottom;
            }
            else
            {
                arrowPrice = Bars.HighPrices[index];
                arrowText = "▼";
                arrowColor = ArrowSellColor;
                vAlign = VerticalAlignment.Top;
            }

            // Draw text as arrow.
            Chart.DrawText(arrowName, arrowText, index, arrowPrice, arrowColor);
            
            // Adjust size.
            var textObj = Chart.FindObject(arrowName) as ChartText;
            if (textObj != null)
            {
                textObj.FontSize = 8 + ArrowSize * 2;
                textObj.VerticalAlignment = vAlign;
                textObj.HorizontalAlignment = HorizontalAlignment.Center;
            }
        }

        private void RemoveCurrentArrow(int index)
        {
            DateTime arrowTime = Bars.OpenTimes[index];
            string arrowName = $"{IndicatorName}-ARROW-{arrowTime:yyyyMMddHHmmss}";
            Chart.RemoveObject(arrowName);
        }

        private void CheckAlerts(int index)
        {
            // Check time constraints.
            if (TriggerCandle == CandleToCheck.ClosedCandle)
            {
                if (MTFBars.OpenTimes[MTFBars.Count - 1] <= LastNotificationTime)
                    return;
            }
            else
            {
                if ((Server.Time - LastNotificationTime).TotalSeconds < WaitTimeNotify)
                    return;
            }
            TradeSignal signal = IsSignal(index);

            if (signal == LastSignal || signal == TradeSignal.Neutral)
                return;

            LastSignal = signal;

            string timeframeStr = MTFBars != Bars ? InputTimeFrame.ToString() : Chart.TimeFrame.ToString();
            string alertMessage = $"{IndicatorName} {Symbol.Name} @ {timeframeStr}: {signal}";

            if (EnableNativeAlerts)
            {
                Notifications.ShowPopup("Coppock Alert", alertMessage, PopupNotificationState.Information);
                Print(alertMessage);
            }

            if (EnableEmailAlerts)
            {
                string subject = $"{IndicatorName} {Symbol.Name} Notification ({timeframeStr}) {signal}";
                string body = $"{Account.BrokerName} - {Account.Number}\n" +
                             $"{IndicatorName} Notification for {Symbol.Name} @ {timeframeStr}\n" +
                             $"Signal: {signal}\n" +
                             $"Time: {Server.Time}";
                
                Notifications.SendEmail(EmailAddress, EmailAddress, subject, body);
            }

            if (EnableSoundAlerts)
            {
                Notifications.PlaySound(SoundType);
            }

            LastNotificationTime = TriggerCandle == CandleToCheck.ClosedCandle ? MTFBars.OpenTimes[MTFBars.Count - 1] : Server.Time;
        }
    }
}