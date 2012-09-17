package es.usal.voronto.view;

import java.awt.BorderLayout;
import java.awt.FlowLayout;

import javax.swing.ComboBoxModel;
import javax.swing.DefaultComboBoxModel;
import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.border.EmptyBorder;
import javax.swing.JLabel;
import javax.swing.JComboBox;
import javax.swing.JList;
import javax.swing.border.LineBorder;

import es.usal.voronto.model.expression.ExpressionData;

import java.awt.Color;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.awt.GridLayout;
import javax.swing.JTextField;

public class DifExpSearch extends JDialog {

	private final JPanel contentPanel = new JPanel();
	private JScrollPane list2;
	private VoronoiVisualization v;
	private JComboBox comboBox;
	private JList group1;
	private JList group2;
	private JTextField textField;

	
	/**
	 * Create the dialog.
	 */
	public DifExpSearch(VoronoiVisualization v) {
		setBounds(100, 100, 412, 461);
		this.v=v;
		getContentPane().setLayout(new BorderLayout());
		contentPanel.setBackground(Color.WHITE);
		contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
		
		getContentPane().add(contentPanel, BorderLayout.CENTER);
		contentPanel.setLayout(null);
		{
			JLabel lblForConditions = new JLabel("terms");
			lblForConditions.setBounds(297, 33, 89, 16);
			contentPanel.add(lblForConditions);
		}
		{
			JLabel lblSearchFor = new JLabel("Search for ");
			lblSearchFor.setBounds(65, 33, 67, 16);
			contentPanel.add(lblSearchFor);
		}
		
		ComboBoxModel regulationModel = 
				new DefaultComboBoxModel(
						new String[] {"up regulated", "down regulated" });
		comboBox = new JComboBox(regulationModel);
		comboBox.setBounds(133, 29, 165, 27);
			
		
		contentPanel.add(comboBox);
		
		JLabel lblOnConditions = new JLabel("on conditions");
		lblOnConditions.setBounds(38, 74, 165, 16);
		contentPanel.add(lblOnConditions);
		
		JLabel lblRespectToConditions = new JLabel("respect to conditions");
		lblRespectToConditions.setBounds(38, 206, 173, 16);
		contentPanel.add(lblRespectToConditions);
		
		list2 = new JScrollPane(); 
		list2.setBounds(188, 206, 196, 118);
		contentPanel.add(list2);
		
		DefaultComboBoxModel model2 = new DefaultComboBoxModel();
		for(String c:v.expData.conditionNames)		model2.addElement(c);
		
		group2 = new JList();
		group2.setModel(model2);
		
		list2.setViewportView(group2);

		
		JScrollPane list1 = new JScrollPane();
		list1.setBounds(188, 74, 198, 120);
		contentPanel.add(list1);
		
		DefaultComboBoxModel model1 = new DefaultComboBoxModel();
		for(String c:v.expData.conditionNames)		model1.addElement(c);
		
		group1 = new JList();
		group1.setModel(model1);
		
		list1.setViewportView(group1);
		
		JLabel lblThreshold = new JLabel("Threshold: ");
		lblThreshold.setBounds(38, 349, 111, 16);
		//lblThreshold.setToolTipText("By default, the standard deviation");
		contentPanel.add(lblThreshold);
		
		textField = new JTextField();
		textField.setBounds(188, 343, 67, 28);
		contentPanel.add(textField);
		textField.setColumns(10);
		textField.setText(1.0+"");

		JPanel buttonPane = new JPanel();
		buttonPane.setBackground(Color.WHITE);
		getContentPane().add(buttonPane, BorderLayout.SOUTH);
		
		JButton okButton = new JButton("OK");
		okButton.setActionCommand("OK");
		okButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				deSearch();
				setVisible(false);
				}
			});
		
		JButton cancelButton = new JButton("Cancel");
		cancelButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				setVisible(false);
				}
			});
		buttonPane.setLayout(new FlowLayout(FlowLayout.CENTER, 5, 5));
		cancelButton.setActionCommand("Cancel");
		buttonPane.add(cancelButton);
		
		buttonPane.add(okButton);
		getRootPane().setDefaultButton(okButton);
	}
	
	private void deSearch()
		{
		v.deSearch(group1.getSelectedIndices(), group2.getSelectedIndices(), comboBox.getSelectedIndex(), new Float(textField.getText()).floatValue());
		}
}
