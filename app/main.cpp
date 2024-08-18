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
#include <QJsonObject>
#include <QJsonDocument>
#include <QRegularExpression>
#include <QFileInfo>

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
        pairTreeWidget->setHeaderLabels(QStringList() << "Pair" << "Forward" << "Reverse");
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
    QString outputFilePath;
    QJsonObject pairs;

    void selectFiles()
    {
        QStringList files = QFileDialog::getOpenFileNames(this, "Select .ab1 Files", "", "AB1 Files (*.ab1)");
        if (files.isEmpty())
            return;

        fileListWidget->clear();
        fileListWidget->addItems(files);
        fileCountLabel->setText(QString("%1 files selected").arg(files.count()));

        updatePairs(files);
    }

    void updatePairs(const QStringList &files)
    {
        pairs = QJsonObject();
        QRegularExpression rx("(.*?)(\\d+)([FR]).*");

        for (const QString &file : files)
        {
            QFileInfo fileInfo(file);
            QString fileName = fileInfo.fileName();
            QRegularExpressionMatch match = rx.match(fileName);
            if (match.hasMatch())
            {
                QString prefix = match.captured(1);
                QString number = match.captured(2);
                QString direction = match.captured(3);
                QString pairKey = prefix + number;
                if (!pairs.contains(pairKey))
                {
                    pairs[pairKey] = QJsonObject({{"F", ""}, {"R", ""}});
                }
                QJsonObject pairObj = pairs[pairKey].toObject();
                pairObj[direction] = file;  // Keep the full path for the file
                pairs[pairKey] = pairObj;
            }
        }

        qDebug() << "Updated pairs:" << QJsonDocument(pairs).toJson(QJsonDocument::Compact);
        updatePairTree();
    }

    void updatePairTree()
    {
        pairTreeWidget->clear();

        for (auto it = pairs.begin(); it != pairs.end(); ++it)
        {
            QJsonObject pairObj = it.value().toObject();
            QTreeWidgetItem *item = new QTreeWidgetItem(pairTreeWidget);
            item->setText(0, it.key());
            item->setText(1, QFileInfo(pairObj["F"].toString()).fileName());
            item->setText(2, QFileInfo(pairObj["R"].toString()).fileName());
        }
    }

    void createConsensus()
    {
        if (pairs.isEmpty())
        {
            QMessageBox::warning(this, "No Pairs Found", "Please select .ab1 files with matching F and R pairs.");
            return;
        }

        QJsonDocument doc(pairs);
        QString pairsJson = doc.toJson(QJsonDocument::Compact);
        qDebug() << "Pairs JSON:" << pairsJson;

        outputFilePath = QDir::tempPath() + "/consensus_output.fasta";
        qDebug() << "Output file path:" << outputFilePath;

        QProcess process;
        process.start("python3", QStringList() << "sequence_aligner.py" << pairsJson << outputFilePath);
        process.waitForFinished(-1);

        QString output(process.readAllStandardOutput());
        qDebug() << "Python output:" << output.trimmed();

        QString error(process.readAllStandardError());
        if (!error.isEmpty()) {
            qDebug() << "Python error output:" << error.trimmed();
        }

        QFile outputFile(outputFilePath);
        if (outputFile.exists()) {
            qDebug() << "Output file exists";
            if (outputFile.size() > 0) {
                qDebug() << "Output file is not empty";
                downloadButton->setEnabled(true);
            } else {
                qDebug() << "Output file is empty";
                QMessageBox::warning(this, "Consensus Creation Failed", "The output file is empty. Please check the Python script for errors.");
            }
        } else {
            qDebug() << "Output file does not exist";
            QMessageBox::warning(this, "Consensus Creation Failed", "The output file was not created. Please check the Python script for errors.");
        }
    }

    void downloadConsensus()
    {
        QString fileName = QFileDialog::getSaveFileName(this, tr("Save Consensus Sequences"), "", tr("FASTA Files (*.fasta)"));
        if (fileName.isEmpty())
            return;

        // Ensure the file has a .fasta extension
        if (!fileName.endsWith(".fasta", Qt::CaseInsensitive))
        {
            fileName += ".fasta";
        }

        if (QFile::copy(outputFilePath, fileName))
        {
            QMessageBox::information(this, "Download Successful", "Consensus sequences saved successfully.");
        }
        else
        {
            QMessageBox::warning(this, "Download Failed", "Failed to save consensus sequences.");
        }
    }
};

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    return a.exec();
}
