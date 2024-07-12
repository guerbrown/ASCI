#include <QApplication>
#include <QMainWindow>
#include <QPushButton>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QProcess>
#include <QDebug>
#include <QListWidget>
#include <QTreeWidget>
#include <QFileDialog>
#include <QLabel>
#include <QSplitter>
#include <QFile>
#include <QTextStream>
#include <QMessageBox>

class MainWindow : public QMainWindow
{
public:
    MainWindow(QWidget *parent = nullptr) : QMainWindow(parent)
    {
        QWidget *centralWidget = new QWidget(this);
        QVBoxLayout *mainLayout = new QVBoxLayout(centralWidget);

        QSplitter *splitter = new QSplitter(Qt::Horizontal);
        mainLayout->addWidget(splitter);

        // Left side
        QWidget *leftWidget = new QWidget;
        QVBoxLayout *leftLayout = new QVBoxLayout(leftWidget);
        
        QHBoxLayout *fileSelectionLayout = new QHBoxLayout();
        QPushButton *selectFilesButton = new QPushButton("Select .ab1 Files", this);
        fileSelectionLayout->addWidget(selectFilesButton);
        
        fileCountLabel = new QLabel("No files selected", this);
        fileSelectionLayout->addWidget(fileCountLabel);
        
        leftLayout->addLayout(fileSelectionLayout);

        fileListWidget = new QListWidget(this);
        leftLayout->addWidget(fileListWidget);

        splitter->addWidget(leftWidget);

        // Right side
        QWidget *rightWidget = new QWidget;
        QVBoxLayout *rightLayout = new QVBoxLayout(rightWidget);
        
        pairTreeWidget = new QTreeWidget(this);
        pairTreeWidget->setHeaderLabels(QStringList() << "Pair" << "Alignment Output");
        rightLayout->addWidget(pairTreeWidget);

        splitter->addWidget(rightWidget);

        createConsensusButton = new QPushButton("Create Sanger Consensus", this);
        mainLayout->addWidget(createConsensusButton);

        downloadButton = new QPushButton("Download Consensus", this);
        downloadButton->setEnabled(false);
        mainLayout->addWidget(downloadButton);

        setCentralWidget(centralWidget);

        connect(selectFilesButton, &QPushButton::clicked, this, &MainWindow::selectFiles);
        connect(createConsensusButton, &QPushButton::clicked, this, &MainWindow::createConsensus);
        connect(downloadButton, &QPushButton::clicked, this, &MainWindow::downloadConsensus);
    }

private:
    QListWidget *fileListWidget;
    QTreeWidget *pairTreeWidget;
    QPushButton *createConsensusButton;
    QPushButton *downloadButton;
    QLabel *fileCountLabel;

    void selectFiles()
    {
        QStringList files = QFileDialog::getOpenFileNames(this, "Select .ab1 Files", "", "AB1 Files (*.ab1)");
        if (files.isEmpty())
            return;

        fileListWidget->clear();
        fileListWidget->addItems(files);
        fileCountLabel->setText(QString("%1 files selected").arg(files.count()));
    }

	   void createConsensus()
	{
	    QStringList inputFiles;
	    for (int i = 0; i < fileListWidget->count(); ++i)
	    {
	        inputFiles << fileListWidget->item(i)->text();
	    }
	
	    if (inputFiles.isEmpty())
	    {
	        QMessageBox::warning(this, "No Files Selected", "Please select .ab1 files first.");
	        return;
	    }
	
	    QProcess process;
	    process.start("python3", QStringList() << "main.py" << inputFiles);
	    process.waitForFinished(-1);
	
	    QString output(process.readAllStandardOutput());
	    qDebug() << "Python output:" << output.trimmed();
	
	    updatePairTree(inputFiles);
	
	    downloadButton->setEnabled(true);
	}

	void updatePairTree(const QStringList &files)
	{
	    pairTreeWidget->clear();
	    QMap<QString, QStringList> pairs;
	
	    QRegExp rx(".*?(\\d+)[FR].*");
	    for (const QString &file : files)
	    {
	        if (rx.exactMatch(file))
	        {
	            QString number = rx.cap(1);
	            pairs[number].append(file);
	        }
	    }
	
	    for (auto it = pairs.begin(); it != pairs.end(); ++it)
	    {
	        QTreeWidgetItem *item = new QTreeWidgetItem(pairTreeWidget);
	        item->setText(0, "Pair " + it.key());
	        item->setExpanded(false);
	
	        QTreeWidgetItem *childItem = new QTreeWidgetItem(item);
	        childItem->setText(0, "Consensus");
	        childItem->setText(1, "Click to load");
	
	        connect(pairTreeWidget, &QTreeWidget::itemClicked, this, [this, it](QTreeWidgetItem *item) {
	            if (item->parent() && item->text(0) == "Consensus")
	            {
	                QString consensusOutput = loadConsensusOutput(it.key());
	                item->setText(1, consensusOutput);
	            }
	        });
	    }
	}

QString loadConsensusOutput(const QString &pairNumber)
{
    // This function should load the consensus output for the given pair number
    // For now, we'll just return a placeholder
    return "Consensus sequence for pair " + pairNumber;
}
    QString loadAlignmentOutput(const QStringList &pairFiles)
    {
        // This is a placeholder. You should implement the actual loading of the alignment output.
        // For now, we'll just return the filenames.
        return pairFiles.join("\n");
    }

    void downloadConsensus()
    {
        QString saveDir = QFileDialog::getExistingDirectory(this, "Select Directory to Save Consensus");
        if (saveDir.isEmpty())
            return;

        // Here you would implement the logic to copy the consensus files to the selected directory
        // For example:
        // QFile::copy("path/to/consensus.fasta", saveDir + "/consensus.fasta");

        qDebug() << "Consensus saved to:" << saveDir;
    }
};

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    return a.exec();
}
