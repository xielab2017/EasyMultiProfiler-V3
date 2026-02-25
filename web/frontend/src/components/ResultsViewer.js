import React, { useState, useEffect } from 'react';
import { 
  Card, Spin, Result, Button, Tabs, Table, Image, 
  Space, Typography, Alert, Progress, Tag, Empty, Row, Col
} from 'antd';
import { 
  ReloadOutlined, DownloadOutlined, FilePdfOutlined, 
  FileExcelOutlined, CheckCircleOutlined, RedoOutlined 
} from '@ant-design/icons';

const { Title, Text } = Typography;
const { TabPane } = Tabs;

function ResultsViewer({ results, taskId, onNewAnalysis }) {
  const [loading, setLoading] = useState(!results);
  const [progress, setProgress] = useState(0);
  const [status, setStatus] = useState('running');

  useEffect(() => {
    if (!results && taskId) {
      const interval = setInterval(async () => {
        try {
          const response = await fetch(`/api/status/${taskId}`);
          const data = await response.json();
          setProgress(data.progress || 0);
          setStatus(data.status);
          
          if (data.status === 'completed') {
            clearInterval(interval);
            setLoading(false);
          } else if (data.status === 'failed') {
            clearInterval(interval);
            setLoading(false);
            setStatus('failed');
          }
        } catch (error) {
          console.error('Error fetching status:', error);
        }
      }, 2000);

      return () => clearInterval(interval);
    }
  }, [results, taskId]);

  if (loading) {
    return (
      <Card style={{ textAlign: 'center', padding: 48 }}>
        <Spin size="large" />
        <Title level={4} style={{ marginTop: 24 }}>分析进行中...⏳</Title>
        <Progress percent={progress} status="active" style={{ maxWidth: 400, margin: '24px auto' }} />
        
        <Text type="secondary">
          {status === 'running' ? '正在执行分析，请稍候...' : '正在初始化...'}
        </Text>
      </Card>
    );
  }

  if (status === 'failed') {
    return (
      <Result
        status="error"
        title="分析失败"
        subTitle="请检查输入数据格式是否正确，或联系管理员"
        extra={[
          <Button type="primary" key="retry" onClick={onNewAnalysis}>
            重新分析
          </Button>
        ]}
      />
    );
  }

  if (!results) {
    return (
      <Empty
        description="暂无分析结果"
        image={Empty.PRESENTED_IMAGE_SIMPLE}
      >
        <Button type="primary" onClick={onNewAnalysis}>开始新分析</Button>
      </Empty>
    );
  }

  const renderPlots = () => {
    if (!results.plots || results.plots.length === 0) {
      return <Empty description="无可用图表" />;
    }

    return (
      <Row gutter={[16, 16]}>
        {results.plots.map((plot, idx) => (
          <Col span={plot.fullWidth ? 24 : 12} key={idx}>
            <Card 
              title={plot.title} 
              size="small"
              extra={
                plot.download_url && (
                  <Button 
                    type="link" 
                    icon={<DownloadOutlined />}
                    href={plot.download_url}
                  >
                    下载
                  </Button>
                )
              }
            >
              <Image 
                src={plot.url} 
                alt={plot.title}
                style={{ width: '100%' }}
                placeholder={<div style={{ height: 300, background: '#f5f5f5' }} />}
              />
              {plot.description && (
                <Text type="secondary" style={{ display: 'block', marginTop: 8 }}>
                  {plot.description}
                </Text>
              )}
            </Card>
          </Col>
        ))}
      </Row>
    );
  };

  const renderTables = () => {
    if (!results.tables || results.tables.length === 0) {
      return <Empty description="无可用表格" />;
    }

    return (
      <Space direction="vertical" style={{ width: '100%' }} size="large">
        {results.tables.map((table, idx) => (
          <Card 
            key={idx}
            title={table.title}
            size="small"
            extra={
              table.download_url && (
                <Button 
                  type="link" 
                  icon={<FileExcelOutlined />}
                  href={table.download_url}
                >
                  下载 CSV
                </Button>
              )
            }
          >
            <Table 
              dataSource={table.data.map((row, i) => ({ ...row, key: i }))}
              columns={table.columns.map(col => ({
                title: col,
                dataIndex: col,
                key: col,
                ellipsis: true
              }))}
              pagination={{ pageSize: 10 }}
              scroll={{ x: 'max-content' }}
              size="small"
            />
          </Card>
        ))}
      </Space>
    );
  };

  const renderStats = () => {
    if (!results.stats) {
      return <Empty description="无统计数据" />;
    }

    const stats = results.stats;
    
    return (
      <Row gutter={[16, 16]}>
        {Object.entries(stats).map(([key, value]) => (
          <Col span={8} key={key}>
            <Card size="small">
              <Text type="secondary">{key}</Text>
              <br />
              <Text strong style={{ fontSize: 24, color: '#1890ff' }}>
                {typeof value === 'number' ? value.toLocaleString() : value}
              </Text>
            </Card>
          </Col>
        ))}
      </Row>
    );
  };

  return (
    <div>
      <Result
        status="success"
        title="分析完成！🎉"
        subTitle={`任务ID: ${taskId}`}
        extra={[
          results.report_url && (
            <Button 
              key="report" 
              type="primary" 
              icon={<FilePdfOutlined />}
              href={results.report_url}
              target="_blank"
            >
              下载完整报告
            </Button>
          ),
          <Button key="new" onClick={onNewAnalysis} icon={<RedoOutlined />}>
            新分析
          </Button>
        ]}
      />

      <Tabs defaultActiveKey="plots" style={{ marginTop: 24 }}>
        <TabPane tab="📊 图表" key="plots">
          {renderPlots()}
        </TabPane>
        
        <TabPane tab="📋 数据表" key="tables">
          {renderTables()}
        </TabPane>
        
        <TabPane tab="📈 统计摘要" key="stats">
          {renderStats()}
        </TabPane>
        
        {results.log && (
          <TabPane tab="📝 分析日志" key="log">
            <pre style={{ 
              background: '#f6ffed', 
              padding: 16, 
              borderRadius: 4,
              maxHeight: 400,
              overflow: 'auto'
            }}>
              {results.log}
            </pre>
          </TabPane>
        )}
      </Tabs>
    </div>
  );
}

export default ResultsViewer;
