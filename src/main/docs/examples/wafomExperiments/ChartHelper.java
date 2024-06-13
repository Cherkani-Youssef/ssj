package wafomExperiments;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import javax.swing.JFrame;
import java.awt.Color;
import java.awt.Font;
import org.jfree.chart.ChartUtilities;
import java.io.File;
import java.io.IOException;



// helper class for plots
public class ChartHelper {

    private String chartTitle;
    private String xAxisLabel;
    private String yAxisLabel;
    private XYSeriesCollection dataset;

    public ChartHelper(String chartTitle, String xAxisLabel, String yAxisLabel, XYSeriesCollection dataset) {
        this.chartTitle = chartTitle;
        this.xAxisLabel = xAxisLabel;
        this.yAxisLabel = yAxisLabel;
        this.dataset = dataset;
    }
    
    
 // New constructor to handle 2D array input
    public ChartHelper(String chartTitle, String xAxisLabel, String yAxisLabel, double[][] data, String[] seriesNames, int k) {
        this.chartTitle = chartTitle;
        this.xAxisLabel = xAxisLabel;
        this.yAxisLabel = yAxisLabel;
        this.dataset = createDatasetFrom2DArray(data, seriesNames, k);
    }

    private XYSeriesCollection createDatasetFrom2DArray(double[][] data, String[] seriesNames , int k) {
        XYSeriesCollection dataset = new XYSeriesCollection();
        for (int i = 0; i < data.length; i++) {
            XYSeries series = new XYSeries(seriesNames != null && seriesNames.length > i ? seriesNames[i] : "Series " + (i + 1));
            for (int j = 0; j < data[i].length; j++) {
                series.add(j +k, data[i][j]);
            }
            dataset.addSeries(series);
        }
        return dataset;
    }

    public void displayChart() {
        JFreeChart chart = createChart();
        styleChart(chart);
        configureDomainAxis(chart);

        ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        JFrame frame = new JFrame(chartTitle);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(chartPanel);
        frame.pack();
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    }

    private JFreeChart createChart() {
        return ChartFactory.createXYLineChart(
            chartTitle,
            xAxisLabel,
            yAxisLabel,
            dataset,
            PlotOrientation.VERTICAL,
            true, // Legend
            true, // Tooltips
            false // URLs
        );
    }

    private void styleChart(JFreeChart chart) {
        XYPlot plot = chart.getXYPlot();
        plot.setBackgroundPaint(Color.white);
//        plot.setDomainGridlinesVisible(true);
        plot.setRangeGridlinesVisible(true);
        plot.setDomainGridlinePaint(new Color(220, 220, 220));
        plot.setRangeGridlinePaint(new Color(220, 220, 220));

        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
        renderer.setSeriesPaint(0, Color.blue.darker()); // NiedXing
		renderer.setSeriesPaint(1, Color.red.darker()); // Sobol
		renderer.setSeriesPaint(2, Color.BLACK); // NiedXing NUS
		renderer.setSeriesPaint(3, Color.yellow.darker()); // Sobol NUS
		renderer.setSeriesPaint(4, Color.cyan); // Explicit NUS
		renderer.setSeriesPaint(5, Color.green.darker()); // Fixed Explicit NUS
//        renderer.setSeriesPaint(0, new Color(31, 119, 180)); // Example color for series 1
//        renderer.setSeriesPaint(1, new Color(255, 127, 14)); // Example color for series 2
        plot.setRenderer(renderer);

        plot.getDomainAxis().setLabelFont(new Font("Tahoma", Font.PLAIN, 12));
        plot.getRangeAxis().setLabelFont(new Font("Tahoma", Font.PLAIN, 12));
        chart.getLegend().setItemFont(new Font("Tahoma", Font.PLAIN, 12));
    }

    private void configureDomainAxis(JFreeChart chart) {
        XYPlot plot = chart.getXYPlot();
        NumberAxis domainAxis = (NumberAxis) plot.getDomainAxis();
        domainAxis.setTickUnit(new NumberTickUnit(1));
        domainAxis.setNumberFormatOverride(new java.text.NumberFormat() {
            private static final long serialVersionUID = 1L;

			@Override
            public StringBuffer format(double number, StringBuffer toAppendTo, java.text.FieldPosition pos) {
                int exp = (int) number;
                return toAppendTo.append("2^").append(exp);
            }

            @Override
            public StringBuffer format(long number, StringBuffer toAppendTo, java.text.FieldPosition pos) {
                return format((double) number, toAppendTo, pos);
            }

            @Override
            public Number parse(String source, java.text.ParsePosition parsePosition) {
                return null; // No parsing is necessary
            }
        });
    }
    
    
    
    
    public void displayScatterPlot() {
        JFreeChart chart = ChartFactory.createScatterPlot(
            chartTitle,
            xAxisLabel,
            yAxisLabel,
            dataset,
            PlotOrientation.VERTICAL,
            true,  // Legend
            true,  // Tooltips
            false  // URLs
        );

        styleScatterPlot(chart);

        ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        JFrame frame = new JFrame(chartTitle);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(chartPanel);
        frame.pack();
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    }

    private void styleScatterPlot(JFreeChart chart) {
        XYPlot plot = chart.getXYPlot();
        plot.setBackgroundPaint(Color.white);
        plot.setDomainGridlinesVisible(true);
        plot.setRangeGridlinesVisible(true);
        plot.setDomainGridlinePaint(Color.lightGray);
        plot.setRangeGridlinePaint(Color.lightGray);

        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
        renderer.setBaseShapesVisible(true);
        renderer.setBaseShapesFilled(true);
        renderer.setDrawSeriesLineAsPath(false);

        // Configure colors or keep it dynamic
        for (int i = 0; i < dataset.getSeriesCount(); i++) {
            renderer.setSeriesPaint(i, Color.getHSBColor((float) i / dataset.getSeriesCount(), 0.85f, 0.85f));
        }
        
        plot.setRenderer(renderer);
    }

    
    public void saveChartAsPNG(String filePath) {
        JFreeChart chart = createChart();
        styleChart(chart);
        configureDomainAxis(chart);

        File outputFile = new File(filePath);
        try {
            // width = 800, height = 600, outputFile = file path to save
            ChartUtilities.saveChartAsPNG(outputFile, chart, 800, 600);
        } catch (IOException e) {
            System.err.println("Problem occurred creating chart.");
            e.printStackTrace();
        }
    }
}
