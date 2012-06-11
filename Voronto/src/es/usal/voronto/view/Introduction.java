package es.usal.voronto.view;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import javax.swing.ComboBoxModel;
import javax.swing.DefaultComboBoxModel;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;

import javax.swing.JLabel;
import javax.swing.filechooser.FileFilter;

import es.usal.voronto.control.Voronto;


/**
* This code was edited or generated using CloudGarden's Jigloo
* SWT/Swing GUI Builder, which is free for non-commercial
* use. If Jigloo is being used commercially (ie, by a corporation,
* company or business for any purpose whatever) then you
* should purchase a license for each developer using Jigloo.
* Please visit www.cloudgarden.com for details.
* Use of Jigloo implies acceptance of these licensing terms.
* A COMMERCIAL LICENSE HAS NOT BEEN PURCHASED FOR
* THIS MACHINE, SO JIGLOO OR THIS CODE CANNOT BE USED
* LEGALLY FOR ANY CORPORATE OR COMMERCIAL PURPOSE.
*/
public class Introduction extends javax.swing.JPanel implements MouseListener{
	/**
	 * 
	 */
	private static final long serialVersionUID = 2715114188900530489L;
	private JLabel jLabel1;
	private JLabel jLabel2;
	private JLabel jLabel3;
	private JLabel jLabel4;
	private JButton jButton1;
	private JLabel jLabel5;
	private JButton jButton2;
	private JLabel jLabel6;
	
	public JComboBox jComboBox1;
	public File expressionFile;
	private File ontologyFile;
	private Voronto parent;
	private JLabel jLabel8;
	private JLabel jLabel9;

	/**
	* Auto-generated main method to display this 
	* JPanel inside a new JFrame.
	*/
		
	public Introduction(Voronto p) {
		super();
		parent=p;
		initGUI();
	}
	
	private void initGUI() {
		try {
			this.setPreferredSize(new java.awt.Dimension(507, 300));
			this.setLayout(null);
			this.setBackground(new java.awt.Color(255,255,255));
			this.setSize(700, 400);
			{
				jLabel1 = new JLabel();
				this.add(jLabel1);
				jLabel1.setText("Mapping gene expression to ontologies with voronoi tessellations");
				jLabel1.setBounds(194, 72, 357, 17);
				jLabel1.setFont(new java.awt.Font("Dialog",2,9));
			}
			{
				jLabel2 = new JLabel();
				this.add(jLabel2);
				jLabel2.setIcon(new ImageIcon(getClass().getClassLoader().getResource("es/usal/voronto/img/voronto.png")));
				jLabel2.setBounds(283, 17, 143, 55);
			}
			{
				jLabel3 = new JLabel();
				this.add(jLabel3);
				jLabel3.setText("1) Select gene expression data");
				jLabel3.setBounds(124, 124, 241, 13);
			}
			{
				jLabel4 = new JLabel();
				this.add(jLabel4);
				jLabel4.setText("Click here for sample data and format specs");
				jLabel4.setBounds(136, 146, 267, 16);
				jLabel4.setFont(new java.awt.Font("Dialog",0,10));
				jLabel4.addMouseListener(this);
				
			}
			{
				jButton1 = new JButton();
				this.add(jButton1);
				jButton1.setText("Select...");
				jButton1.setBounds(390, 124, 100, 24);
				jButton1.setBackground(new java.awt.Color(255,255,255));
				jButton1.setSize(98, 22);
				jButton1.addActionListener(new java.awt.event.ActionListener()	{
					public void actionPerformed(java.awt.event.ActionEvent e) 
						{
						JFileChooser selecFile = new JFileChooser();
						ExpressionDataFilter edf=new ExpressionDataFilter();
						selecFile.addChoosableFileFilter(edf);
						selecFile.setFileFilter(edf);
						//selecFile.setCurrentDirectory(new File("data"));
						
						String path=null;
						try{
							BufferedReader br=new BufferedReader(new FileReader(parent.expressionPath));
							path=br.readLine();
							}catch(Exception ex){System.out.println("pathFile non existing yet");}
						if(path!=null)				
							selecFile.setCurrentDirectory(new File(path));
						else						
							selecFile.setCurrentDirectory(new File("."));
						
						int returnval = selecFile.showDialog((Component)e.getSource(), "Load expression data");
						
						if(returnval == JFileChooser.APPROVE_OPTION) 
							{
							expressionFile = selecFile.getSelectedFile();
							jLabel8.setText(expressionFile.getAbsolutePath());
							jButton2.setEnabled(true);
							}
					}
					});
			}
			
			{
				jLabel5 = new JLabel();
				this.add(jLabel5);
				jLabel5.setText("2) Select ontology");
				jLabel5.setBounds(124, 218, 162, 14);
			}
			{
				ComboBoxModel jComboBox1Model = 
						new DefaultComboBoxModel(
								new String[] { "KEGG", "GO slim BP","GO slim CC", "GO full BP", "REACTOME", "Custom..." });
				jComboBox1 = new JComboBox();
				this.add(jComboBox1);
				jComboBox1.setModel(jComboBox1Model);
				jComboBox1.setBounds(380, 214, 130, 22);
				jComboBox1.setBackground(new java.awt.Color(255,255,255));
				jComboBox1.setSelectedIndex(0);
				jComboBox1.addActionListener(new java.awt.event.ActionListener() {

					//private File ontologyFile;

					
					@Override
					public void actionPerformed(ActionEvent e) 
						{
						if(jComboBox1.getSelectedIndex()==5)
							{
							JFileChooser selecFile = new JFileChooser();
							XMLDataFilter edf=new XMLDataFilter();
							selecFile.addChoosableFileFilter(edf);
							selecFile.setFileFilter(edf);
							
							String path=null;
							try{
								BufferedReader br=new BufferedReader(new FileReader(parent.ontologyPath));
								path=br.readLine();
								}catch(Exception ex){System.out.println("pathFile non existing yet");}
							if(path!=null)				
								selecFile.setCurrentDirectory(new File(path));
							else						
								selecFile.setCurrentDirectory(new File("."));
							
							int returnval = selecFile.showDialog((Component)e.getSource(), "Load ontology");
							
							if(returnval == JFileChooser.APPROVE_OPTION) 
								{
								ontologyFile = selecFile.getSelectedFile();
								jLabel9.setText(ontologyFile.getAbsolutePath());
								}
							}
						else
							{
							ontologyFile = null;
							jLabel9.setText("");
							}
						}
					});
			}
			{
				jLabel6 = new JLabel();
				this.add(jLabel6);
				jLabel6.setText("3) Launch application!");
				jLabel6.setBounds(124, 263, 162, 19);
				
			}
			{
				jButton2 = new JButton();
				this.add(jButton2);
				jButton2.setText("Launch");
				jButton2.setBackground(new java.awt.Color(255,255,255));
				jButton2.setBounds(392, 265, 100, 22);
				jButton2.setEnabled(false);
				jButton2.addActionListener(new java.awt.event.ActionListener() {
					public void actionPerformed(java.awt.event.ActionEvent e) 
						{
						if(expressionFile==null)	
							{
							if(ontologyFile==null)	parent.launch(null, null, jComboBox1.getSelectedIndex());
							else					parent.launch(null, ontologyFile.getPath(), jComboBox1.getSelectedIndex());
							}
						else
							{
							if(ontologyFile==null)	parent.launch(expressionFile.getPath(), null, jComboBox1.getSelectedIndex());
							else					parent.launch(expressionFile.getPath(), ontologyFile.getPath(), jComboBox1.getSelectedIndex());
							}
						}
					});
			}
			{
				jLabel8 = new JLabel();
				this.add(jLabel8);
				jLabel8.setText("No file selected");
				jLabel8.setBounds(136, 166, 400, 15);
				jLabel8.setFont(new java.awt.Font("Dialog",0,10));
				jLabel8.setForeground(new Color(150,150,255));
			}
			{
				jLabel9 = new JLabel();
				this.add(jLabel9);
				jLabel9.setText("");
				jLabel9.setBounds(136, 238, 397, 15);
				jLabel9.setFont(new java.awt.Font("Dialog",0,10));
				jLabel9.setForeground(new Color(150,150,255));
			}

		} catch (Exception e) {
			e.printStackTrace();
		}
	}


/**
 * Expression file format filter for JFileChooser
 * @author Rodrigo	Santamaria
 */
private static class ExpressionDataFilter extends FileFilter{
    static final String getExtension(File f) {
        String ext = null;
        String s = f.getName();
        int i = s.lastIndexOf('.');

        if (i > 0 &&  i < s.length() - 1) {
            ext = s.substring(i+1).toLowerCase();
        }
        return ext;
    }
    
  
   /* static final String getExtension(String s) 
    	{
        String ext = null;
        int i = s.lastIndexOf('.');

        if (i > 0 &&  i < s.length() - 1)          ext = s.substring(i+1).toLowerCase();
        return ext;
    	}*/
    
    /**
     * Decides if a file is acceptable as Microarray file
     * @param	f	File to check
     * @return	true if the file extension is txt, false otherwise
     */
    public boolean accept(File f) 
    	{
        if (f.isDirectory()) 							       return true;

        String extension = getExtension(f);
        if (extension != null) 
        	{
            if (extension.equals("txt"))                  return true;
            else								          return false;
            }

        return false;
    	}

    //The description of this filter
    /**
     * Returns the description of Microarray files
     * @return	A brief String description of expected files for Microarray.
     */
    public String getDescription() {
        return "Tab delimited data (.txt)";
    }
}


/**
 * Expression file format filter for JFileChooser
 * @author Rodrigo	Santamaria
 */
private static class XMLDataFilter extends FileFilter{
    static final String getExtension(File f) {
        String ext = null;
        String s = f.getName();
        int i = s.lastIndexOf('.');

        if (i > 0 &&  i < s.length() - 1) {
            ext = s.substring(i+1).toLowerCase();
        }
        return ext;
    }
    
  /*
    static final String getExtension(String s) 
    	{
        String ext = null;
        int i = s.lastIndexOf('.');

        if (i > 0 &&  i < s.length() - 1)          ext = s.substring(i+1).toLowerCase();
        return ext;
    	}*/
    
    /**
     * Decides if a file is acceptable as Microarray file
     * @param	f	File to check
     * @return	true if the file extension is txt, false otherwise
     */
    public boolean accept(File f) 
    	{
        if (f.isDirectory()) 							       return true;

        String extension = getExtension(f);
        if (extension != null) 
        	{
            if (extension.equals("xml"))                  return true;
            else								          return false;
            }

        return false;
    	}

    //The description of this filter
    /**
     * Returns the description of Microarray files
     * @return	A brief String description of expected files for Microarray.
     */
    public String getDescription() {
        return "XML ontology data (.xml)";
    }
}

public void reset()
	{
	//jLabel8.setText("");
	jLabel9.setText("");
	if(this.expressionFile==null)	jButton2.setEnabled(false);
	}

@Override
public void mouseClicked(MouseEvent arg0) {
	// TODO Auto-generated method stub
	
}

@Override
public void mouseEntered(MouseEvent arg0) {
	// TODO Auto-generated method stub
	
}

@Override
public void mouseExited(MouseEvent arg0) {
	// TODO Auto-generated method stub
	
}

@Override
public void mousePressed(MouseEvent arg0) {
	// TODO Auto-generated method stub
	
}

@Override
public void mouseReleased(MouseEvent arg0) {
	if(arg0.getSource().getClass()==JLabel.class)
		{
		if((JLabel)(arg0.getSource())==jLabel4)
			{
			try{
				java.awt.Desktop.getDesktop().browse(java.net.URI.create("http://vis.usal.es/~visusal/voronto/voronto/Help.html"));
			}catch(IOException ioe){ioe.printStackTrace();}
			}
		}
}
}